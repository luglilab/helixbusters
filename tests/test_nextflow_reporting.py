"""Exercise production DSL2 wiring with stub outputs, without genomic inputs."""

import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).parents[1]
STUBS = {
    "CHECK_ENVIRONMENT": "echo '{}' > environment.json",
    "EXTRACT_UMI": "touch ${meta.sample}.umi.fastq.gz",
    "MAP_READS": """touch ${meta.sample}.q${params.mapq}.bam ${meta.sample}.q${params.mapq}.bam.bai
    touch ${meta.sample}.all.bam ${meta.sample}.all.bam.bai ${meta.sample}.mapping.json
    echo ${task.cpus} > ${meta.sample}.cpu.log""",
    "DEDUPLICATE": """touch ${meta.sample}.families.tsv ${meta.sample}.counts.bed ${meta.sample}.dedup.json
    touch ${meta.sample}.molecules.bed ${meta.sample}.sites.tsv""",
    "SAMPLE_QC": """touch ${meta.sample}.summary.json ${meta.sample}.chrom.sizes.json ${meta.sample}.qc.tsv
    touch sample__${meta.sample}.filtered.flagstat.txt sample__${meta.sample}.all.flagstat.txt
    touch ${meta.sample}.ends.raw.bw ${meta.sample}.ends.CPM.bw ${meta.sample}.coverage.CPM.bw
    touch ${meta.sample}.ends.provenance.json ${meta.sample}.ends.counts.bed""",
    "CONDITION_QC": """touch ${group}.condition.summary.json ${group}.ends.counts.bed
    touch ${group}.ends.raw.bw ${group}.ends.CPM.bw ${group}.ends.mean.CPM.bw ${group}.ends.provenance.json
    touch ${group}.condition.filtered.bam ${group}.condition.filtered.bam.bai
    touch condition__${group}.filtered.flagstat.txt
    echo ${summaries.size()} > ${group}.condition.qc.tsv""",
    "MULTIQC": """touch multiqc_report.html reporting_versions.txt
    mkdir multiqc_data
    echo ${samples.size()} > helixbusters_samples.tsv
    echo ${conditions.size()} > helixbusters_conditions.tsv""",
}


@unittest.skipUnless(shutil.which("nextflow"), "Nextflow is not on PATH")
class TestNextflowReporting(unittest.TestCase):
    def test_three_samples_two_conditions_and_integer_cpu_requests(self):
        with tempfile.TemporaryDirectory(prefix="helix-nextflow-") as directory:
            root = Path(directory)
            # Only add stubs to the copied workflow. Production scripts remain
            # in place and are compiled; stub data never enter real analyses.
            source = (ROOT / "main.nf").read_text()
            pattern = re.compile(r"process (\w+) \{(.*?)(?=\nprocess |\nworkflow \{)", re.S)

            def add_stub(match):
                name, block = match.groups()
                script = block.index("    script:")
                return f"process {name} {{" + block[:script] + '    stub:\n    """\n    ' + STUBS[name] + '\n    """\n\n' + block[script:]

            (root / "main.nf").write_text(pattern.sub(add_stub, source))
            shutil.copytree(ROOT / "conf", root / "conf")
            shutil.copy2(ROOT / "nextflow.config", root / "nextflow.config")
            (root / "stub.config").write_text(
                "process {\n memory = '256 MB'\n"
                "withLabel: small { memory = '256 MB' }\n"
                "withLabel: mapping { memory = '256 MB' }\n"
                "withLabel: reporting { memory = '256 MB' }\n}\n"
                "executor.cpus = 16\nexecutor.queueSize = 1\n"
                "timeline.enabled = false\nreport.enabled = false\ntrace.enabled = false\ndag.enabled = false\n")
            (root / "input.fastq").write_text("@test\nACGT\n+\nIIII\n")
            (root / "reference.json").write_text("{}\n")
            (root / "samples.tsv").write_text(
                "sample\tgroup\treplicate\tbarcode\tfastq\n" + "".join(
                    f"{sample}\t{group}\t{replicate}\tACGT\t{root / 'input.fastq'}\n"
                    for sample, group, replicate in (("a", "treated", "1"), ("b", "treated", "2"), ("c", "control", "1"))))
            result = subprocess.run([
                "nextflow", "run", "main.nf", "-stub-run", "-profile", "local", "-c", "stub.config",
                "--manifest", "samples.tsv", "--genome", "hg38", "--reference_config", "reference.json",
                "--map_threads", "5", "--sort_threads", "5", "--coverage_bin_size", "50",
                "--outdir", "outputs"], cwd=root, env={**os.environ, "NXF_OFFLINE": "true"},
                capture_output=True, text=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            outputs = root / "outputs"
            for sample in ("a", "b", "c"):
                self.assertEqual((outputs / "SingleReplicate" / sample / "mapping" / f"{sample}.cpu.log").read_text().strip(), "11")
                self.assertTrue((outputs / "SingleReplicate" / sample / "bigwig" / f"{sample}.ends.CPM.bw").is_file())
            for group, n in (("treated", 2), ("control", 1)):
                self.assertEqual((outputs / "MergedReplicate" / group / "qc" / f"{group}.condition.qc.tsv").read_text().strip(), str(n))
                self.assertTrue((outputs / "MergedReplicate" / group / "mapping" / f"{group}.condition.filtered.bam").is_file())
            self.assertEqual((outputs / "MultiQC" / "helixbusters_samples.tsv").read_text().strip(), "3")
            self.assertEqual((outputs / "MultiQC" / "helixbusters_conditions.tsv").read_text().strip(), "2")


if __name__ == "__main__":
    unittest.main()
