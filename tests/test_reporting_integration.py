"""Real BAM -> library tracks -> condition merge -> MultiQC on a toy genome."""

import importlib.util
import gzip
import json
import math
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import pysam

from helixbusters.deduplication import deduplicate_bam
from helixbusters.mapping import filter_alignments, map_sample
from tests.bam_helpers import read, write_bam


ROOT = Path(__file__).parents[1]
AVAILABLE = all(shutil.which(name) for name in ("samtools", "bamCoverage", "multiqc")) and importlib.util.find_spec("pyBigWig")


@unittest.skipUnless(AVAILABLE, "Reporting integration requires samtools, deepTools, MultiQC and pyBigWig")
class TestReportingIntegration(unittest.TestCase):
    @unittest.skipUnless(shutil.which("bwa"), "BWA is not on PATH")
    def test_real_bwa_alignment_to_library_report_and_bigwig(self):
        import pyBigWig
        with tempfile.TemporaryDirectory(prefix="helix-bwa-report-") as directory:
            root = Path(directory)
            rng = random.Random(51)
            sequence = "".join(rng.choice("ACGT") for _ in range(4000))
            fasta = root / "reference.fa"
            fasta.write_text(">chrSynthetic\n" + sequence + "\n")
            subprocess.run(["bwa", "index", str(fasta)], check=True, capture_output=True)
            reads = root / "reads.fastq.gz"
            with gzip.open(reads, "wt") as handle:
                for i, umi in enumerate(("AAAAAAAA", "AAAAAAAA", "CCCCCCCC")):
                    seq = sequence[300:385]
                    handle.write(f"@read{i}_{umi}\n{seq}\n+\n{'I' * len(seq)}\n")
            output = root / "output"
            mapped = map_sample("library", reads, fasta, output, aligner="bwa", threads=1, sort_threads=0)
            count_path, dedup_path = output / "library.counts.bed", output / "library.dedup.json"
            qc = deduplicate_bam(mapped["BamFilteredPath"], output / "library.families.tsv", count_path,
                                 method="directional", umi_length=8, min_mapq=20, output_qc=dedup_path)
            self.assertEqual(qc["deduplicated_molecules"], 2)
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/post_mapping.py"), "sample",
                "--sample", "library", "--group", "condition", "--replicate", "1",
                "--mapping", mapped["MappingQC"], "--dedup", str(dedup_path),
                "--all-bam", mapped["BamAllPath"], "--filtered-bam", mapped["BamFilteredPath"],
                "--counts", str(count_path), "--threads", "1"], cwd=output,
                capture_output=True, text=True, timeout=90)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            bw = pyBigWig.open(str(output / "library.ends.raw.bw"))
            try:
                self.assertEqual(bw.values("chrSynthetic", 300, 301), [2.0])
            finally:
                bw.close()

    def test_independent_deduplication_real_tracks_and_merged_qc(self):
        import pyBigWig

        with tempfile.TemporaryDirectory(prefix="helix-reporting-") as directory:
            root = Path(directory)
            summaries, counts, headers, bams = [], [], [], []
            fixtures = {
                "a": [read(umi="AAAAAAAA"), read(umi="AAAAAAAA"), read(umi="CCCCCCCC"),
                      read(start=200, umi="GGGGGGGG", flag=16), read(start=300, mapq=10),
                      read(ref=-1, start=-1, flag=4, cigar=None)],
                "b": [read(umi="AAAAAAAA"), read(start=200, umi="CCCCCCCC", flag=16),
                      read(start=300, mapq=10), read(ref=-1, start=-1, flag=4, cigar=None)],
            }

            def run(*arguments, cwd):
                result = subprocess.run(arguments, cwd=cwd, capture_output=True, text=True, timeout=90)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            for replicate, (sample, reads) in enumerate(fixtures.items(), 1):
                out = root / sample
                out.mkdir()
                all_bam = write_bam(out / f"{sample}.all.bam", reads, header={
                    "HD": {"VN": "1.6", "SO": "coordinate"},
                    "SQ": [{"SN": "chrSynthetic", "LN": 1000}]})
                pysam.index(str(all_bam))
                filtered = out / f"{sample}.filtered.bam"
                mapping = filter_alignments(all_bam, filtered, 20)
                pysam.index(str(filtered))
                mapping.update(sample=sample, parameters={"aligner": "synthetic BAM fixture", "min_mapq": 20})
                mapping_path = out / f"{sample}.mapping.json"
                mapping_path.write_text(json.dumps(mapping))
                count_path, dedup_path = out / f"{sample}.counts.bed", out / f"{sample}.dedup.json"
                deduplicate_bam(filtered, out / f"{sample}.families.tsv", count_path,
                                method="directional", umi_length=8, min_mapq=20, output_qc=dedup_path)
                run(sys.executable, str(ROOT / "scripts/post_mapping.py"), "sample",
                    "--sample", sample, "--group", "condition", "--replicate", str(replicate),
                    "--mapping", str(mapping_path), "--dedup", str(dedup_path),
                    "--all-bam", str(all_bam), "--filtered-bam", str(filtered), "--counts", str(count_path),
                    "--threads", "1", "--bin-size", "50", cwd=out)
                summaries.append(out / f"{sample}.summary.json")
                headers.append(out / f"{sample}.chrom.sizes.json")
                counts.append(count_path)
                bams.append(filtered)
                bw = pyBigWig.open(str(out / f"{sample}.coverage.CPM.bw"))
                try:
                    self.assertTrue(any(math.isfinite(value) and value > 0
                                        for value in bw.values("chrSynthetic", 100, 150)))
                finally:
                    bw.close()

            group = root / "condition"
            group.mkdir()
            run(sys.executable, str(ROOT / "scripts/post_mapping.py"), "group",
                "--group", "condition", "--summaries", *map(str, summaries),
                "--headers", *map(str, headers), "--counts", *map(str, counts),
                "--bams", *map(str, bams), "--threads", "1", cwd=group)
            with pysam.AlignmentFile(str(group / "condition.condition.filtered.bam"), "rb") as merged:
                self.assertEqual(sum(1 for _ in merged), 6)
            for suffix, expected in (("raw", 3), ("CPM", 600000), ("mean.CPM", 583333.333333)):
                bw = pyBigWig.open(str(group / f"condition.ends.{suffix}.bw"))
                try:
                    self.assertAlmostEqual(bw.values("chrSynthetic", 100, 101)[0], expected, delta=0.1)
                finally:
                    bw.close()
            condition_summary = group / "condition.condition.summary.json"
            pooled = json.loads(condition_summary.read_text())["metrics"]
            self.assertEqual(pooled["deduplicated_molecules"], 5)
            self.assertEqual(pooled["excluded_low_mapq"], 2)
            report = root / "report"
            report.mkdir()
            run(sys.executable, str(ROOT / "scripts/post_mapping.py"), "report",
                "--samples", *map(str, summaries), "--conditions", str(condition_summary), cwd=report)
            for directory in (root / "a", root / "b", group):
                for path in directory.glob("*.txt"):
                    shutil.copy2(path, report / path.name)
            run("multiqc", ".", "--filename", "multiqc_report.html", "--outdir", ".",
                "--data-dir", "--cl-config", "data_dir_name: multiqc_data", cwd=report)
            self.assertTrue((report / "multiqc_data").is_dir())
            self.assertTrue(any((report / "multiqc_data").iterdir()))
            html = (report / "multiqc_report.html").read_text()
            for text in ("Helixbusters samples", "Helixbusters conditions", "Samtools"):
                self.assertIn(text, html)


if __name__ == "__main__":
    unittest.main()
