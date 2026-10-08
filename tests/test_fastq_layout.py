import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/inspect_fastq_layout.py"
SPEC = importlib.util.spec_from_file_location("inspect_fastq_layout", SCRIPT)
layout = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(layout)


class TestFastqLayout(unittest.TestCase):
    def test_positions_orientations_and_bounded_gzip_read(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq.gz"
            barcode = "CTCACACG"
            sequences = ["AACCGGTT" + barcode + "GATTACA",
                         "AA" + layout.reverse_complement(barcode) + "GATTACA",
                         "TTTTTTTT" + barcode + "GATTACA"]
            with gzip.open(path, "wt") as handle:
                for i, seq in enumerate(sequences):
                    handle.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
                handle.write("@malformed_after_limit\n")
            result = layout.inspect_fastq(path, [barcode], 2, 80, 20, 8, 8)
            self.assertEqual(result["reads_examined"], 2)
            hits = {hit["orientation"]: hit for hit in result["barcode_search"]}
            self.assertEqual(hits["forward"]["positions"][0]["offset_0based"], 8)
            self.assertEqual(hits["reverse_complement"]["positions"][0]["offset_0based"], 2)
            self.assertEqual(result["forward_barcode_at_tested_offset"][barcode]["exact_reads"], 1)
            with self.assertRaisesRegex(ValueError, "record 4"):
                list(layout.fastq_sequences(path, 4))

    def test_manifest_resolves_shared_input_without_altering_it(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "samples.tsv"
            manifest.write_text("sample\tgroup\treplicate\tbarcode\tfastq\n"
                                "a\tCHRONIC\tR1\tACGT\treads.fastq.gz\n"
                                "b\tCHRONIC\tR2\tTGCA\treads.fastq.gz\n")
            original = manifest.read_bytes()
            rows = layout.read_manifest(manifest)
            self.assertEqual(rows[0]["fastq"], rows[1]["fastq"])
            self.assertTrue(Path(rows[0]["fastq"]).is_absolute())
            self.assertEqual(manifest.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
