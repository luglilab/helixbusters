import hashlib
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "prepare_igenome", Path(__file__).parents[1] / "scripts/prepare_igenome.py")
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class TestReferencePreparation(unittest.TestCase):
    def test_resume_retains_partial_and_hashes_complete_file(self):
        with tempfile.TemporaryDirectory() as root:
            target = Path(root) / "archive.partial"
            target.write_bytes(b"first")
            commands = []

            def run(command, **kwargs):
                commands.append(command)
                with target.open("ab") as handle:
                    handle.write(b"next")
                return subprocess.CompletedProcess(command, 18 if len(commands) == 1 else 0,
                                                   stderr="connection interrupted")

            with patch.object(prepare.shutil, "which", return_value="curl"), \
                 patch.object(prepare.subprocess, "run", side_effect=run), \
                 patch.object(prepare.time, "sleep"):
                digest = prepare.download("https://example.org/archive", target, "100K", 2)
            self.assertEqual(target.read_bytes(), b"firstnextnext")
            self.assertEqual(digest, hashlib.sha256(target.read_bytes()).hexdigest())
            self.assertIn("--continue-at", commands[0])
            self.assertIn("100K", commands[0])

    def test_http_failure_stops_and_preserves_partial(self):
        with tempfile.TemporaryDirectory() as root:
            target = Path(root) / "archive.partial"
            target.write_bytes(b"saved")
            with patch.object(prepare.shutil, "which", return_value="curl"), \
                 patch.object(prepare.subprocess, "run", return_value=
                              subprocess.CompletedProcess([], 22, stderr="HTTP 404")) as run:
                with self.assertRaisesRegex(ValueError, "HTTP 404"):
                    prepare.download("https://example.org/archive", target)
            self.assertEqual(run.call_count, 1)
            self.assertEqual(target.read_bytes(), b"saved")

    def test_shared_versioned_index_requires_nonempty_files(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            version = directory / "version1"
            version.mkdir()
            prefix = version / "genome.fa"
            for suffix in prepare.INDEX_SUFFIXES:
                Path(str(prefix) + suffix).write_bytes(b"index")
            self.assertEqual(prepare.find_bwa_prefix(directory), prefix)
            Path(str(prefix) + ".sa").write_bytes(b"")
            with self.assertRaisesRegex(ValueError, "nonempty"):
                prepare.find_bwa_prefix(directory)


if __name__ == "__main__":
    unittest.main()
