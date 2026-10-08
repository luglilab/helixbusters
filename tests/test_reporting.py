import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "reporting", Path(__file__).parents[1] / "helixbusters/reporting.py")
reporting = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reporting)


class TestReporting(unittest.TestCase):
    def test_multiqc_hides_noncanonical_contigs_without_changing_idxstats(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            idxstats = root / 'sample.all.idxstats.txt'
            original = ('chr1\t248956422\t100\t0\n1\t248956422\t20\t0\n'
                        'chrM\t16569\t10\t0\nchrUn_GL000220v1\t1000\t5\t0\n'
                        'chr17_KI270729v1_random\t1000\t3\t0\n')
            idxstats.write_text(original)
            with patch.object(reporting.Path, 'glob', return_value=[idxstats]), \
                 patch.object(reporting, 'write_table'), \
                 patch.object(reporting, 'write_json') as write_json:
                reporting.report_content([], [], 'hg38')
            config = write_json.call_args_list[0].args[1]
            self.assertEqual(config['samtools_idxstats_ignore'],
                             ['chr17_KI270729v1_random', 'chrM', 'chrUn_GL000220v1'])
            self.assertEqual(idxstats.read_text(), original)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.header = [("chr2", 100), ("chr1", 100)]

    def bed(self, name, text):
        path = self.root / name
        path.write_text(text)
        return path

    def test_pooled_and_equal_weight_mean_differ_with_library_depth(self):
        first = self.bed("a.bed", "chr2\t1\t2\t2\n")
        second = self.bed("b.bed", "chr2\t1\t2\t1\nchr1\t5\t6\t7\n")
        rows = list(reporting.merge_end_counts([first, second], [2, 8], self.header))
        self.assertEqual(rows[0][:4], ("chr2", 1, 2, 3))
        self.assertAlmostEqual(rows[0][4], 300000)
        self.assertAlmostEqual(rows[0][5], 562500)
        self.assertEqual(rows[1][:4], ("chr1", 5, 6, 7))
        self.assertAlmostEqual(sum(row[4] for row in rows), 1e6)
        self.assertAlmostEqual(sum(row[5] for row in rows), 1e6)

    def test_zero_library_has_undefined_mean_and_no_pseudocount(self):
        first = self.bed("a.bed", "")
        second = self.bed("b.bed", "chr2\t1\t2\t2\n")
        rows = list(reporting.merge_end_counts([first, second], [0, 2], self.header))
        self.assertEqual(rows[0][4], 1e6)
        self.assertIsNone(rows[0][5])
        self.assertEqual(list(reporting.merge_end_counts([first], [0], self.header)), [])

    def test_bad_coordinates_order_or_qc_total_fail(self):
        for text, total in [("chr2\t1\t3\t2\n", 2),
                            ("chrX\t1\t2\t2\n", 2),
                            ("chr2\t100\t101\t2\n", 2),
                            ("chr1\t5\t6\t1\nchr2\t1\t2\t1\n", 2),
                            ("chr2\t1\t2\t2\n", 3)]:
            with self.subTest(text=text):
                bed = self.bed("bad.bed", text)
                with self.assertRaises(ValueError):
                    list(reporting.merge_end_counts([bed], [total], self.header))

    def record(self, sample, replicate, primary, mapped, retained, molecules):
        mapping = {"alignment_records": primary, "primary_records": primary,
                   "mapped_primary_records": mapped, "retained_records": retained,
                   "excluded_by_reason": {"mitochondrial": mapped - retained},
                   "parameters": {"aligner": "bwa", "min_mapq": 20}, "genome_filter": {"genome": "hg38"}}
        dedup = {"accepted_reads": retained, "deduplicated_molecules": molecules,
                 "duplicate_reads": retained - molecules, "parameters": {"method": "directional"}}
        return {"sample": sample, "group": "condition", "replicate": replicate,
                "mapping": mapping, "deduplication": dedup, "metrics": reporting.metrics(mapping, dedup)}

    def test_group_uses_ratio_of_sums_and_preserves_replicate_variation(self):
        records = [self.record("a", "1", 10, 10, 8, 5), self.record("b", "2", 90, 45, 40, 30)]
        summaries, headers = [], []
        for record in records:
            summary = self.root / f'{record["sample"]}.json'
            summary.write_text(json.dumps(record))
            summaries.append(summary)
            header = self.root / f'{record["sample"]}.header.json'
            header.write_text(json.dumps(self.header))
            headers.append(header)
        with patch.object(reporting, "write_json") as write_json, \
             patch.object(reporting, "write_table"), patch.object(reporting, "end_tracks") as tracks:
            reporting.group_qc("condition", summaries, ["a.bed", "b.bed"], headers)
        result = write_json.call_args.args[1]["metrics"]
        self.assertEqual(result["mapping_rate_pct"], 55)
        self.assertEqual(result["replicate_mean_mapping_rate_pct"], 75)
        self.assertEqual(result["deduplicated_molecules"], 35)
        self.assertEqual(tracks.call_args.args[1], [5, 30])
        records[1]["replicate"] = "1"
        summaries[1].write_text(json.dumps(records[1]))
        with self.assertRaisesRegex(ValueError, "biological replicate"):
            reporting.group_qc("condition", summaries, ["a.bed", "b.bed"], headers)

    def test_empty_denominators_and_unsafe_names(self):
        row = self.record("a", "1", 0, 0, 0, 0)["metrics"]
        self.assertIsNone(row["mapping_rate_pct"])
        self.assertIsNone(row["umi_duplicate_pct"])
        for name in ("../x", "a b", "a'", "", "x/y"):
            with self.assertRaises(ValueError):
                reporting.validate_label(name)

    def test_missing_worker_tools_fail_before_reference_access(self):
        with patch.object(reporting.shutil, "which", return_value=None):
            with self.assertRaisesRegex(ValueError, "Missing worker executables.*bamCoverage.*multiqc"):
                reporting.check_environment("absent.json", "hg38", "bwa")

    @unittest.skipUnless(importlib.util.find_spec("pyBigWig"), "pyBigWig is not installed")
    def test_real_bigwig_roundtrip_and_zero_library(self):
        import pyBigWig
        first = self.bed("a.bed", "chr2\t1\t2\t2\n")
        second = self.bed("b.bed", "chr2\t1\t2\t1\nchr1\t5\t6\t7\n")
        prefix = self.root / "pooled"
        reporting.end_tracks([first, second], [2, 8], self.header, prefix, mean=True)
        for suffix, expected in (("raw", 3), ("CPM", 300000), ("mean.CPM", 562500)):
            bw = pyBigWig.open(f"{prefix}.ends.{suffix}.bw")
            try:
                self.assertEqual(bw.chroms(), dict(self.header))
                self.assertAlmostEqual(bw.values("chr2", 1, 2)[0], expected)
            finally:
                bw.close()
        original = Path(f"{prefix}.ends.raw.bw").read_bytes()
        with self.assertRaisesRegex(ValueError, "already exist"):
            reporting.end_tracks([first, second], [2, 8], self.header, prefix, mean=True)
        self.assertEqual(Path(f"{prefix}.ends.raw.bw").read_bytes(), original)
        empty = self.bed("empty.bed", "")
        reporting.end_tracks([empty], [0], self.header, self.root / "empty", mean=True)
        self.assertFalse((self.root / "empty.ends.mean.CPM.bw").exists())
        bw = pyBigWig.open(str(self.root / "empty.ends.CPM.bw"))
        try:
            self.assertEqual(bw.header()["nBasesCovered"], 0)
        finally:
            bw.close()

    @unittest.skipUnless(shutil.which("multiqc"), "MultiQC is not on PATH")
    def test_multiqc_parses_both_custom_sections(self):
        record = self.record("a", "1", 10, 10, 8, 5)
        sample = self.root / "sample.summary.json"
        sample.write_text(json.dumps(record))
        condition = self.root / "condition.summary.json"
        condition.write_text(json.dumps({"group": "condition", "metrics": record["metrics"]}))
        original_json, original_table = reporting.write_json, reporting.write_table
        with patch.object(reporting, "write_json", side_effect=lambda path, data: original_json(self.root / path, data)), \
             patch.object(reporting, "write_table", side_effect=lambda path, data: original_table(self.root / path, data)):
            reporting.report_content([sample], [condition])
        result = subprocess.run(["multiqc", ".", "--filename", "report.html", "--outdir", "."],
                                cwd=self.root, capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        html = (self.root / "report.html").read_text()
        self.assertIn("Helixbusters samples", html)
        self.assertIn("Helixbusters conditions", html)


if __name__ == "__main__":
    unittest.main()
