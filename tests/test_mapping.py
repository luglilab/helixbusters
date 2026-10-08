import gzip
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
import pysam

from helixbusters.core import Helixbusters
from helixbusters.deduplication import deduplicate_bam
from helixbusters.mapping import (MappingError, build_alignment_command, filter_alignments,
                                 map_sample, run_alignment_pipe, validate_index)
from tests.bam_helpers import read, write_bam


class MappingFixtures(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="helix mapping ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def executable(self, name, body):
        path = self.root / name
        path.write_text(f"#!{sys.executable}\n" + body)
        path.chmod(0o755)
        return str(path)

    def index(self, extension="bwa"):
        prefix = self.root / "reference genome"
        suffixes = ((".amb", ".ann", ".bwt", ".pac", ".sa") if extension == "bwa" else
                    tuple(s + extension for s in (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")))
        for suffix in suffixes:
            Path(str(prefix) + suffix).write_text("fixture index; not an actual reference")
        return prefix


class TestMappingCommands(MappingFixtures):
    def test_bwa_single_and_paired_argument_vectors(self):
        command = build_alignment_command("bwa", "bwa", "/reference genome", "/read 1;$.fq.gz",
                                          "sample one", 4)
        self.assertEqual(command[-2:], ["/reference genome", "/read 1;$.fq.gz"])
        self.assertIn("@RG\\tID:sample one\\tSM:sample one", command)
        command = build_alignment_command("bwa", "bwa", "ref", "r1", "sample", 4, "r2")
        self.assertEqual(command[-3:], ["ref", "r1", "r2"])

    def test_bowtie_modes_and_mates(self):
        command = build_alignment_command("bowtie2", "bowtie2", "ref", "r1", "s", 2)
        self.assertIn("--end-to-end", command)
        self.assertEqual(command[-2:], ["-U", "r1"])
        command = build_alignment_command("bowtie2", "bowtie2", "ref", "r1", "s", 2,
                                          "r2", "local")
        self.assertIn("--local", command)
        self.assertEqual(command[-4:], ["-1", "r1", "-2", "r2"])

    def test_complete_index_sets(self):
        prefix = self.index()
        self.assertEqual(validate_index(prefix, "bwa"), str(prefix))
        Path(str(prefix) + ".sa").unlink()
        with self.assertRaises(FileNotFoundError):
            validate_index(prefix, "bwa")
        for extension in (".bt2", ".bt2l"):
            prefix = self.index(extension)
            self.assertEqual(validate_index(prefix, "bowtie2"), str(prefix))
            for path in self.root.glob("*" + extension):
                path.unlink()

    def test_bad_options_fail_before_running_programs(self):
        for options in ({"threads": 0}, {"sort_threads": -1}, {"sort_memory": "2;rm"},
                        {"min_mapq": 255}, {"sample": "../sample"}):
            kwargs = dict(sample="sample", read1="missing", genome_index="missing", output_dir=self.root)
            kwargs.update(options)
            with self.subTest(options=options), self.assertRaises(ValueError):
                map_sample(**kwargs)


class TestAlignmentFiltering(MappingFixtures):
    def test_filter_boundaries_flags_and_qc_accounting(self):
        records = [read(mapq=20), read(mapq=19), read(mapq=255), read(flag=256),
                   read(flag=2048), read(flag=512), read(flag=1024), read(cigar="3S20M"),
                   read(flag=65), read(flag=145), read(start=-1, ref=-1, flag=4, cigar=None)]
        source = write_bam(self.root / "input.bam", records)
        target = self.root / "filtered.bam"
        qc = filter_alignments(source, target, 20)
        self.assertEqual(qc["retained_records"], 5)
        self.assertEqual(qc["retained_ambiguous_five_prime"], 1)
        self.assertEqual(qc["excluded_by_reason"], {"unmapped": 1, "secondary": 1,
                         "supplementary": 1, "qc_fail": 1, "unknown_mapq": 1, "low_mapq": 1})
        self.assertEqual(qc["alignment_records"], qc["retained_records"] + sum(qc["excluded_by_reason"].values()))
        self.assertEqual(qc["primary_records"], 9)
        self.assertEqual(qc["mapped_primary_records"], 8)
        self.assertEqual(qc["primary_mapped_mapq_histogram"]["19"], 1)
        with pysam.AlignmentFile(target, "rb") as bam:
            retained = list(bam)
        self.assertTrue(any(r.is_duplicate for r in retained))
        self.assertEqual(sum(r.is_paired for r in retained), 2)

    def test_never_overwrites_filter_input(self):
        source = write_bam(self.root / "input.bam", [read()])
        original = source.read_bytes()
        with self.assertRaises(ValueError):
            filter_alignments(source, source)
        self.assertEqual(source.read_bytes(), original)


class TestCheckedPipe(MappingFixtures):
    def test_mapper_failure_is_detected_when_consumer_succeeds(self):
        mapper = [sys.executable, "-c", "import sys; print('partial SAM'); sys.exit(7)"]
        sorter = [sys.executable, "-c", "import sys; sys.stdin.read()"]
        with self.assertRaisesRegex(MappingError, "Aligner exited 7"):
            run_alignment_pipe(mapper, sorter, self.root / "map.log", self.root / "sort.log")

    def test_sort_failure_stops_running_producer(self):
        mapper = [sys.executable, "-c", "import time; time.sleep(30)"]
        sorter = [sys.executable, "-c", "import sys; sys.stderr.write('disk full'); sys.exit(9)"]
        with self.assertRaisesRegex(MappingError, "sort exited 9"):
            run_alignment_pipe(mapper, sorter, self.root / "map.log", self.root / "sort.log")
        self.assertEqual((self.root / "sort.log").read_text(), "disk full")


class TestMappingPublication(MappingFixtures):
    def setUp(self):
        super().setUp()
        self.prefix = self.index()
        self.fastq = self.root / "reads with spaces.fastq.gz"
        with gzip.open(self.fastq, "wt") as handle:
            handle.write("@read_AAAAAAAA\n" + "A" * 20 + "\n+\n" + "I" * 20 + "\n")
        fixture_bam = write_bam(self.root / "fixture.bam", [read(start=200), read(start=100)], sort=False)
        with pysam.AlignmentFile(fixture_bam, "rb") as bam:
            sam = str(bam.header) + "".join(r.to_string() + "\n" for r in bam)
        self.mapper = self.executable("mapper", "import sys\n"
            "if len(sys.argv) == 1: print('fixture mapper version 1'); sys.exit(0)\n"
            f"sys.stdout.write({sam!r})\n")
        self.sorter = self.executable("samtools", "import sys, pysam\n"
            "if sys.argv[1] == '--version': print('fixture sorter version 1'); sys.exit(0)\n"
            "pysam.sort(*sys.argv[2:])\n")
        self.options = dict(sample="sample", read1=self.fastq, genome_index=self.prefix,
                            output_dir=self.root / "out", threads=1, sort_threads=0,
                            aligner_executable=self.mapper, samtools_executable=self.sorter)

    def test_real_pipe_sort_filter_and_indexes(self):
        outputs = map_sample(**self.options)
        for key in ("BamAllPath", "BamFilteredPath"):
            self.assertTrue(Path(outputs[key] + ".bai").is_file())
            with pysam.AlignmentFile(outputs[key], "rb") as bam:
                self.assertEqual([r.reference_start for r in bam.fetch()], [100, 200])
        qc = json.loads(Path(outputs["MappingQC"]).read_text())
        self.assertEqual(qc["retained_records"], 2)
        self.assertEqual(qc["inputs"]["read1"], str(self.fastq))
        self.assertEqual(qc["parameters"]["min_mapq"], 20)

    def test_failed_rerun_keeps_previous_bams_and_qc(self):
        outputs = map_sample(**self.options)
        old = {key: Path(outputs[key]).read_bytes() for key in ("BamAllPath", "BamFilteredPath", "MappingQC")}
        self.options["aligner_executable"] = self.executable("broken", "import sys\nsys.exit(7)\n")
        with self.assertRaises(MappingError):
            map_sample(**self.options)
        for key, content in old.items():
            self.assertEqual(Path(outputs[key]).read_bytes(), content)
        self.assertEqual(list(Path(self.options["output_dir"]).glob(".sample.mapping-*")), [])

    def test_index_failure_does_not_publish(self):
        with patch("helixbusters.mapping.pysam.index", side_effect=RuntimeError("index failed")):
            with self.assertRaisesRegex(RuntimeError, "index failed"):
                map_sample(**self.options)
        self.assertFalse((Path(self.options["output_dir"]) / "sample.all.bam").exists())
        self.assertTrue((Path(self.options["output_dir"]) / "sample.bwa.log").exists())


class TestMappingIntegrationWithCore(MappingFixtures):
    def test_explicit_trimmed_paths_and_metadata_are_used(self):
        reads = self.root / "custom.fastq.gz"
        reads.touch()
        hb = Helixbusters("unused", "human", 1, "index")
        hb.modality = "single-end"
        hb.infofile = pd.DataFrame([{"Sample": "s", "OutputPath": str(self.root),
                                    "PathReadForwardTrimmed": str(reads)}])
        with patch("helixbusters.core.map_sample", return_value={"MappingQC": "qc.json"}) as call:
            hb.run_bowtie2_mapping(quality=30, threads=3, bowtie2_mode="local")
        self.assertEqual(call.call_args.args[1], str(reads))
        self.assertEqual(call.call_args.kwargs["min_mapq"], 30)
        self.assertEqual(call.call_args.kwargs["bowtie2_mode"], "local")
        self.assertEqual(hb.infofile.loc[0, "MappingQC"], "qc.json")

    def test_missing_later_sample_is_caught_before_mapping_first(self):
        reads = self.root / "exists.fastq.gz"
        reads.touch()
        hb = Helixbusters("unused", "human", 1, "index")
        hb.modality = "single-end"
        hb.infofile = pd.DataFrame([{"Sample": name, "OutputPath": str(self.root),
                                    "PathReadForwardTrimmed": str(path)}
                                   for name, path in (("one", reads), ("two", self.root / "missing"))])
        with patch("helixbusters.core.map_sample") as call:
            with self.assertRaises(FileNotFoundError):
                hb.run_bwa_mapping()
            call.assert_not_called()


@unittest.skipUnless(shutil.which("bowtie2") and shutil.which("bowtie2-build") and shutil.which("samtools"),
                     "Optional real Bowtie2/samtools test: add executables to PATH")
class TestRealBowtie2(MappingFixtures):
    def test_gzipped_reads_mapping_to_deduplication(self):
        rng = random.Random(51)
        sequence = "".join(rng.choice("ACGT") for _ in range(4000))
        fasta = self.root / "reference.fa"
        fasta.write_text(">chrSynthetic\n" + sequence + "\n")
        prefix = str(self.root / "reference_index")
        subprocess.run(["bowtie2-build", "--quiet", str(fasta), prefix], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        forward = sequence[300:385]
        reverse = sequence[900:985].translate(str.maketrans("ACGT", "TGCA"))[::-1]
        reads = self.root / "trimmed.fastq.gz"
        with gzip.open(reads, "wt") as handle:
            for i, (seq, umi) in enumerate([(forward, "AAAAAAAA"), (forward, "AAAAAAAA"),
                                           (forward, "CCCCCCCC"), (reverse, "GGGGGGGG")]):
                handle.write(f"@read{i}_{umi}\n{seq}\n+\n{'I' * len(seq)}\n")
        outputs = map_sample("synthetic", reads, prefix, self.root / "out", aligner="bowtie2",
                             threads=1, sort_threads=0)
        qc = json.loads(Path(outputs["MappingQC"]).read_text())
        self.assertEqual(qc["retained_records"], 4)
        molecules = self.root / "molecules.bed"
        counts = deduplicate_bam(outputs["BamFilteredPath"], self.root / "families.tsv",
                                 self.root / "counts.bed", output_molecules=molecules)
        self.assertEqual(counts["deduplicated_molecules"], 3)
        rows = [line.split("\t") for line in molecules.read_text().splitlines()]
        self.assertEqual([(int(row[1]), row[5]) for row in rows], [(300, "+"), (300, "+"), (984, "-")])


if __name__ == "__main__":
    unittest.main()
