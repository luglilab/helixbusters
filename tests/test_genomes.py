import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
import pysam

from helixbusters.core import Helixbusters
from helixbusters.deduplication import deduplicate_bam
from helixbusters.genomes import (GenomeFilter, canonical_lengths, genome_species,
                                  load_reference_config, make_genome_filter, normalize_genome)
from helixbusters.mapping import filter_alignments, map_sample
from tests.bam_helpers import read, write_bam


def genome_header(genome, prefix=True):
    lengths = canonical_lengths(genome)
    return {"HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": name if prefix else name[3:], "LN": length}
                   for name, length in lengths.items()] +
                  [{"SN": "chrM" if prefix else "MT", "LN": 16569},
                   {"SN": "chrUn_scaffold", "LN": 10000}]}


class GenomeFixtures(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def blacklist(self, text="chr1\t100\t200\n", compressed=False):
        path = self.root / ("blacklist.bed.gz" if compressed else "blacklist.bed")
        if compressed:
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
        return path


class TestGenomeSelection(GenomeFixtures):
    def test_all_aliases_and_canonical_sets(self):
        for name, alias, species, count in [("hg19", "GRCh37", "human", 24),
                                            ("hg38", "GRCh38", "human", 24),
                                            ("mm10", "GRCm38", "mouse", 21),
                                            ("mm39", "GRCm39", "mouse", 21)]:
            self.assertEqual(normalize_genome(alias), name)
            self.assertEqual(genome_species(alias), species)
            self.assertEqual(len(canonical_lengths(alias)), count)
            self.assertNotIn("chrM", canonical_lengths(alias))
            self.assertIn("chrX", canonical_lengths(alias))
            self.assertIn("chrY", canonical_lengths(alias))
        with self.assertRaises(ValueError):
            normalize_genome("mm9")

    def test_select_build_requires_a_matching_blacklist(self):
        path = self.blacklist()
        with self.assertRaisesRegex(ValueError, "does not match"):
            GenomeFilter("mm39", path, "mm10")
        with self.assertRaisesRegex(ValueError, "explicitly"):
            GenomeFilter("hg38", path, None)
        with self.assertRaisesRegex(ValueError, "required"):
            GenomeFilter("hg38", None, "hg38")
        with self.assertRaises(FileNotFoundError):
            GenomeFilter("hg38", self.root / "missing.bed", "hg38")
        self.assertIsNone(make_genome_filter())
        with self.assertRaises(ValueError):
            make_genome_filter(blacklist_bed=path)

    def test_species_is_inferred_and_conflicts_fail(self):
        hb = Helixbusters("samples.csv", genome="GRCm39", genome_index="ref",
                          blacklist_bed=self.blacklist(), blacklist_genome="mm39")
        self.assertEqual(hb.species, "mouse")
        self.assertEqual(hb.genome, "mm39")
        with self.assertRaisesRegex(ValueError, "Species does not match"):
            Helixbusters("samples.csv", species="human", genome="mm39")

    def test_genome_header_length_and_name_checks(self):
        filters = GenomeFilter("hg38", self.blacklist(), "hg38")
        for prefix in (True, False):
            header = genome_header("hg38", prefix)
            filters.validate_header([s["SN"] for s in header["SQ"]], [s["LN"] for s in header["SQ"]])
        header = genome_header("hg19")
        with self.assertRaisesRegex(ValueError, "genome mismatch"):
            filters.validate_header([s["SN"] for s in header["SQ"]], [s["LN"] for s in header["SQ"]])
        with self.assertRaisesRegex(ValueError, "Missing canonical"):
            filters.validate_header(["NC_000001.11"], [248956422])
        with self.assertRaisesRegex(ValueError, "multiple aliases"):
            filters.validate_header(["chr1", "1"], [248956422, 248956422])

    def test_wrong_mouse_build_is_rejected(self):
        filters = GenomeFilter("mm39", self.blacklist(), "mm39")
        header = genome_header("mm10")
        with self.assertRaisesRegex(ValueError, "genome mismatch"):
            filters.validate_header([s["SN"] for s in header["SQ"]], [s["LN"] for s in header["SQ"]])


class TestBlacklist(GenomeFixtures):
    def test_half_open_boundaries_and_aliases(self):
        filters = GenomeFilter("hg38", self.blacklist(), "hg38")
        self.assertFalse(filters.overlaps("chr1", 80, 100))
        self.assertTrue(filters.overlaps("1", 99, 101))
        self.assertTrue(filters.overlaps("chr1", 199, 201))
        self.assertFalse(filters.overlaps("chr1", 200, 220))
        self.assertFalse(filters.overlaps("chr2", 100, 200))

    def test_merging_compression_comments_and_provenance(self):
        path = self.blacklist("# comment\ntrack name=test\nchr1\t150\t250\n"
                              "1\t100\t200\nchr1\t250\t300\nchrUn_scaffold\t1\t20\n", True)
        filters = GenomeFilter("GRCh38", path, "hg38")
        meta = filters.metadata()
        self.assertEqual(meta["blacklist_input_intervals"], 4)
        self.assertEqual(meta["blacklist_merged_intervals"], 1)
        self.assertEqual(meta["blacklist_noncanonical_intervals_ignored"], 1)
        self.assertEqual(meta["blacklist_sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
        self.assertTrue(filters.overlaps("chr1", 299, 301))

    def test_invalid_beds_fail_instead_of_silently_disabling_filter(self):
        for text in ("", "chr1\t-1\t100\n", "chr1\t20\t20\n", "chr1\t30\t20\n",
                     "chr1\tstart\tend\n", "chr1\t0\t999999999\n",
                     "chrUn_scaffold\t0\t100\n", "NC_000001.11\t0\t100\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                GenomeFilter("hg38", self.blacklist(text), "hg38")


class TestGenomeFiltering(GenomeFixtures):
    def test_all_four_builds_filter_and_feed_deduplication(self):
        for genome in ("hg19", "hg38", "mm10", "mm39"):
            with self.subTest(genome=genome):
                header = genome_header(genome)
                canonical = len(canonical_lengths(genome))
                records = [read(start=80), read(start=90), read(start=100), read(start=200),
                           read(start=200, flag=16), read(start=300, ref=canonical - 2),
                           read(start=300, ref=canonical - 1), read(ref=canonical), read(ref=canonical + 1)]
                source = write_bam(self.root / "input.bam", records, header=header)
                filters = GenomeFilter(genome, self.blacklist(), genome)
                output = self.root / "filtered.bam"
                qc = filter_alignments(source, output, genome_filter=filters)
                self.assertEqual(qc["retained_records"], 5)
                self.assertEqual(qc["excluded_by_reason"], {"blacklist": 2, "mitochondrial": 1, "noncanonical": 1})
                self.assertEqual(qc["genome_filter"]["genome"], genome)
                stats = deduplicate_bam(output, self.root / "families.tsv", self.root / "counts.bed")
                self.assertEqual(stats["deduplicated_molecules"], 5)
                self.assertNotIn("chrM", (self.root / "counts.bed").read_text())

    def test_ensembl_names_are_preserved(self):
        header = genome_header("hg38", prefix=False)
        source = write_bam(self.root / "input.bam", [read(start=90), read(start=200), read(ref=24)], header=header)
        output = self.root / "filtered.bam"
        qc = filter_alignments(source, output, genome_filter=GenomeFilter("hg38", self.blacklist(), "hg38"))
        self.assertEqual(qc["excluded_by_reason"], {"blacklist": 1, "mitochondrial": 1})
        with pysam.AlignmentFile(output) as bam:
            self.assertEqual([r.reference_name for r in bam], ["1"])

    def test_mitochondrial_aliases_and_mouse_chr20(self):
        for alias in ("chrM", "chrMT", "MT", "M"):
            header = genome_header("mm39")
            header["SQ"][-2]["SN"] = alias
            header["SQ"].append({"SN": "chr20", "LN": 10000})
            source = write_bam(self.root / "input.bam", [read(ref=21), read(ref=23)], header=header)
            qc = filter_alignments(source, self.root / "out.bam",
                                   genome_filter=GenomeFilter("mm39", self.blacklist(), "mm39"))
            self.assertEqual(qc["excluded_by_reason"], {"mitochondrial": 1, "noncanonical": 1})

    def test_gaps_and_clips_are_not_aligned_bases(self):
        header = genome_header("hg38")
        # Blacklist is [100,200). A D/N spanning the entire region is not aligned
        # read sequence; an M inside the interval is. The 5' end alone is not used.
        records = [read(start=80, cigar="20M100D20M"), read(start=80, cigar="20M100N20M"),
                   read(start=200, cigar="10S20M"), read(start=80, cigar="121M")]
        source = write_bam(self.root / "input.bam", records, header=header)
        qc = filter_alignments(source, self.root / "out.bam",
                               genome_filter=GenomeFilter("hg38", self.blacklist(), "hg38"))
        self.assertEqual(qc["retained_records"], 3)
        self.assertEqual(qc["excluded_by_reason"], {"blacklist": 1})

    def test_wrong_header_does_not_touch_output(self):
        output = self.root / "out.bam"
        output.write_bytes(b"previous output")
        source = write_bam(self.root / "input.bam", [read()], header=genome_header("hg19"))
        with self.assertRaisesRegex(ValueError, "genome mismatch"):
            filter_alignments(source, output, genome_filter=GenomeFilter("hg38", self.blacklist(), "hg38"))
        self.assertEqual(output.read_bytes(), b"previous output")

    def test_missing_blacklist_fails_before_starting_mapper(self):
        with patch("helixbusters.mapping.run_alignment_pipe") as run:
            with self.assertRaises(FileNotFoundError):
                map_sample("s", "reads", "index", self.root, genome="mm39",
                           blacklist_bed=self.root / "missing.bed", blacklist_genome="mm39")
            run.assert_not_called()


class TestReferenceCatalog(GenomeFixtures):
    def config(self):
        bed = self.blacklist()
        path = self.root / "references.json"
        path.write_text(json.dumps({"schema_version": 1, "references": {
            genome: {"bwa_index": genome + "/bwa/genome", "bowtie2_index": genome + "/bt2/genome",
                     "blacklist": {"genome": genome, "path": bed.name}}
            for genome in ("hg19", "hg38", "mm10", "mm39")}}))
        return path

    def test_selection_resolves_paths_and_aligner(self):
        config = self.config()
        for genome, alias in (("hg19", "GRCh37"), ("hg38", "GRCh38"), ("mm10", "GRCm38"), ("mm39", "GRCm39")):
            hb = Helixbusters.from_reference_config("samples.csv", alias, config, aligner="bowtie2")
            self.assertEqual(hb.genome, genome)
            self.assertEqual(hb.genome_index, str(self.root / genome / "bt2/genome"))
            self.assertEqual(hb.blacklist_bed, str(self.root / "blacklist.bed"))
            with self.assertRaisesRegex(ValueError, "selected a bowtie2 index"):
                hb.run_bwa_mapping()

    def test_core_passes_filters_to_mapping(self):
        hb = Helixbusters.from_reference_config("samples.csv", "hg38", self.config())
        reads = self.root / "reads.fastq.gz"
        reads.touch()
        hb.modality = "single-end"
        hb.infofile = pd.DataFrame([{"Sample": "s", "OutputPath": str(self.root),
                                    "PathReadForwardTrimmed": str(reads)}])
        with patch("helixbusters.core.map_sample", return_value={}) as run:
            hb.run_bwa_mapping()
        self.assertEqual(run.call_args.kwargs["genome"], "hg38")
        self.assertEqual(run.call_args.kwargs["blacklist_genome"], "hg38")
        self.assertEqual(run.call_args.kwargs["blacklist_bed"], str(self.root / "blacklist.bed"))

    def test_mismatched_catalog_blacklist_fails(self):
        config = self.config()
        data = json.loads(config.read_text())
        data["references"]["mm39"]["blacklist"]["genome"] = "mm10"
        config.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "does not match"):
            Helixbusters.from_reference_config("samples.csv", "mm39", config)


if __name__ == "__main__":
    unittest.main()
