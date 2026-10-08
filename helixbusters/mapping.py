"""Checked mapper execution, coordinate sorting and BLISS alignment QC."""

from collections import Counter
from pathlib import Path
import json
import os
import re
import shutil
import signal
import subprocess
import tempfile

import pysam

from helixbusters.deduplication import five_prime_position
from helixbusters.genomes import make_genome_filter


class MappingError(RuntimeError):
    """A mapper or sorting process failed; diagnostic logs are preserved."""


def validate_sample_name(sample):
    if (not isinstance(sample, str) or not sample.strip() or sample in {".", ".."}
            or any(c in sample for c in "/\\\t\r\n")
            or any(ord(c) < 32 for c in sample)):
        raise ValueError("sample must be a nonempty filename-safe name without control characters")


def validate_index(prefix, aligner):
    """Require a complete classic BWA or Bowtie2 small/large index."""
    prefix = str(Path(prefix).resolve())
    if aligner == "bwa":
        candidates = [[prefix + suffix for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]]
    elif aligner == "bowtie2":
        candidates = [[prefix + suffix + extension
                       for suffix in (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")]
                      for extension in (".bt2", ".bt2l")]
    else:
        raise ValueError("aligner must be 'bwa' or 'bowtie2'")
    for files in candidates:
        if all(Path(path).is_file() and Path(path).stat().st_size > 0 for path in files):
            return prefix
    raise FileNotFoundError(f"Incomplete or missing {aligner} index: {prefix}")


def _executable(name):
    executable = shutil.which(os.fspath(name))
    if executable is None:
        raise FileNotFoundError(f"Required executable not found: {name}")
    return str(Path(executable).absolute())


def build_alignment_command(aligner, executable, index, read1, sample, threads,
                            read2=None, bowtie2_mode="end-to-end"):
    """Build argument vectors; never pass filenames through a shell."""
    if aligner == "bwa":
        rg = f"@RG\\tID:{sample}\\tSM:{sample}"
        return [executable, "mem", "-v", "1", "-t", str(threads), "-R", rg,
                index, read1] + ([read2] if read2 is not None else [])
    if aligner != "bowtie2":
        raise ValueError("aligner must be 'bwa' or 'bowtie2'")
    if bowtie2_mode not in {"end-to-end", "local"}:
        raise ValueError("bowtie2_mode must be 'end-to-end' or 'local'")
    command = [executable, "--" + bowtie2_mode, "-x", index, "-p", str(threads),
               "--rg-id", sample, "--rg", f"SM:{sample}"]
    return command + (["-U", read1] if read2 is None else ["-1", read1, "-2", read2])


def _stop(process):
    if process is not None and process.poll() is None:
        def send(sig):
            try:
                if os.name == "posix":
                    os.killpg(process.pid, sig)
                else:
                    process.send_signal(sig)
            except ProcessLookupError:
                pass

        # Bowtie2 is a wrapper: terminate its group, including the aligner child.
        send(signal.SIGTERM)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            if os.name == "posix":
                send(signal.SIGKILL)
            else:
                process.kill()
            process.wait()


def run_alignment_pipe(mapper_command, sort_command, mapper_log, sort_log):
    """Check both exit statuses and reap children, including on interruption."""
    mapper = sorter = None
    try:
        with open(mapper_log, "wb") as map_errors, open(sort_log, "wb") as sort_errors:
            mapper = subprocess.Popen(mapper_command, stdout=subprocess.PIPE, stderr=map_errors,
                                      start_new_session=(os.name == "posix"))
            try:
                sorter = subprocess.Popen(sort_command, stdin=mapper.stdout,
                                          stdout=subprocess.DEVNULL, stderr=sort_errors,
                                          start_new_session=(os.name == "posix"))
            finally:
                mapper.stdout.close()
            sort_status = sorter.wait()
            if sort_status:
                _stop(mapper)
                raise MappingError(f"samtools sort exited {sort_status}; see {sort_log} and {mapper_log}")
            mapper_status = mapper.wait()
            if mapper_status:
                raise MappingError(f"Aligner exited {mapper_status}; see {mapper_log}")
    finally:
        _stop(sorter)
        _stop(mapper)


def filter_alignments(input_bam, output_bam, min_mapq=20, *, genome_filter=None):
    """Keep primary mapped QC-passing records with known MAPQ >= threshold.

    Count alignment records, not molecules or paired fragments. Keep both mates,
    coordinate-only duplicate flags, and 5-prime clipped records for downstream
    explicit end selection. The complete BAM remains available for re-filtering.
    """
    if type(min_mapq) is not int or not 0 <= min_mapq <= 254:
        raise ValueError("min_mapq must be an integer between 0 and 254")
    if Path(input_bam).resolve() == Path(output_bam).resolve():
        raise ValueError("Input and output BAM paths must be distinct")
    counts = Counter({"alignment_records": 0, "primary_records": 0,
                      "mapped_primary_records": 0, "retained_records": 0,
                      "retained_ambiguous_five_prime": 0})
    excluded, mapq = Counter(), Counter()
    with pysam.AlignmentFile(str(input_bam), "rb") as source:
        if genome_filter is not None:
            genome_filter.validate_header(source.references, source.lengths)
        with pysam.AlignmentFile(str(output_bam), "wb", template=source) as target:
            for read in source.fetch(until_eof=True):
                counts["alignment_records"] += 1
                if not read.is_secondary and not read.is_supplementary:
                    counts["primary_records"] += 1
                    if not read.is_unmapped:
                        counts["mapped_primary_records"] += 1
                        mapq[read.mapping_quality] += 1
                reason = None
                if read.is_unmapped:
                    reason = "unmapped"
                elif read.is_secondary:
                    reason = "secondary"
                elif read.is_supplementary:
                    reason = "supplementary"
                elif read.is_qcfail:
                    reason = "qc_fail"
                elif read.mapping_quality == 255:
                    reason = "unknown_mapq"
                elif read.mapping_quality < min_mapq:
                    reason = "low_mapq"
                elif genome_filter is not None:
                    reason = genome_filter.exclusion_reason(read)
                if reason:
                    excluded[reason] += 1
                else:
                    target.write(read)
                    counts["retained_records"] += 1
                    counts["retained_ambiguous_five_prime"] += five_prime_position(read) is None
    return dict(counts, excluded_by_reason=dict(sorted(excluded.items())),
                primary_mapped_mapq_histogram={str(k): v for k, v in sorted(mapq.items())},
                genome_filter=genome_filter.metadata() if genome_filter else None)


def _version(executable, aligner=False):
    # Classic bwa prints its version in help and returns nonzero without args.
    try:
        result = subprocess.run([executable] if aligner else [executable, "--version"],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, errors="replace", timeout=10)
        return result.stdout.strip()[:2000]
    except (OSError, subprocess.TimeoutExpired) as error:
        return f"Version unavailable: {error}"


def map_sample(sample, read1, genome_index, output_dir, *, read2=None, aligner="bwa",
               min_mapq=20, threads=10, sort_threads=1, sort_memory="768M",
               bowtie2_mode="end-to-end", aligner_executable=None,
               samtools_executable=None, genome=None, blacklist_bed=None,
               blacklist_genome=None):
    """Map one sample, publish verified BAMs/BAIs and return output paths.

    BAMs and report are staged until alignment, filtering and both indexes
    succeed. Diagnostic logs describe the latest attempt, including failures.
    """
    validate_sample_name(sample)
    genome_filter = make_genome_filter(genome, blacklist_bed, blacklist_genome)
    if type(min_mapq) is not int or not 0 <= min_mapq <= 254:
        raise ValueError("min_mapq must be an integer between 0 and 254")
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    if type(sort_threads) is not int or sort_threads < 0:
        raise ValueError("sort_threads must be a nonnegative integer")
    if not isinstance(sort_memory, str) or not re.fullmatch(r"[1-9][0-9]*[KMG]", sort_memory):
        raise ValueError("sort_memory must be an integer with K, M or G suffix")
    reads = [str(Path(path).resolve()) for path in (read1, read2) if path is not None]
    if read1 is None:
        raise ValueError("read1 is required")
    for path in reads:
        if not Path(path).is_file():
            raise FileNotFoundError(f"Trimmed FASTQ file not found: {path}")
        if aligner == "bowtie2" and "," in path:
            raise ValueError("Bowtie2 treats commas as file-list separators; rename the FASTQ path")
    if len(reads) == 2 and reads[0] == reads[1]:
        raise ValueError("Read1 and read2 must be distinct files")
    index = validate_index(genome_index, aligner)
    mapper = _executable(aligner_executable or aligner)
    samtools = _executable(samtools_executable or "samtools")
    command = build_alignment_command(aligner, mapper, index, reads[0], sample, threads,
                                      reads[1] if len(reads) == 2 else None, bowtie2_mode)
    output_dir = Path(output_dir).resolve()
    names = {"BamAllPath": f"{sample}.all.bam", "BamFilteredPath": f"{sample}.q{min_mapq}.bam",
             "MappingQC": f"{sample}.mapping.json", "MappingLog": f"{sample}.{aligner}.log",
             "SortLog": f"{sample}.sort.log"}
    outputs = {key: str(output_dir / name) for key, name in names.items()}
    produced = list(outputs.values()) + [outputs[key] + ".bai" for key in ("BamAllPath", "BamFilteredPath")]
    if set(reads) & {str(Path(path).resolve()) for path in produced}:
        raise ValueError("Output paths must not overwrite input FASTQs")
    if genome_filter is not None and str(genome_filter.path) in produced:
        raise ValueError("Output paths must not overwrite the blacklist")
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{sample}.mapping-", dir=output_dir) as tmp:
        work = Path(tmp)
        bam_all = work / names["BamAllPath"]
        bam_filtered = work / names["BamFilteredPath"]
        sort_command = [samtools, "sort", "-@", str(sort_threads), "-m", sort_memory,
                        "-T", str(work / "sort"), "-o", str(bam_all), "-"]
        run_alignment_pipe(command, sort_command, outputs["MappingLog"], outputs["SortLog"])
        qc = filter_alignments(bam_all, bam_filtered, min_mapq, genome_filter=genome_filter)
        # Index creation validates coordinate ordering and consumes each BAM.
        for path in (bam_all, bam_filtered):
            pysam.index(str(path))
        qc["sample"] = sample
        qc["inputs"] = {"read1": reads[0], "read2": reads[1] if len(reads) == 2 else None,
                        "genome_index": index}
        qc["parameters"] = {"aligner": aligner, "min_mapq": min_mapq, "threads": threads,
                            "sort_threads": sort_threads, "sort_memory_per_thread": sort_memory,
                            "alignment_mode": bowtie2_mode if aligner == "bowtie2" else "bwa-mem"}
        qc["commands"] = {"aligner": command, "sort": sort_command}
        qc["versions"] = {"aligner": _version(mapper, aligner == "bwa"),
                          "samtools": _version(samtools), "pysam": pysam.__version__}
        (work / names["MappingQC"]).write_text(json.dumps(qc, indent=2) + "\n")
        for name in (names["BamAllPath"], names["BamAllPath"] + ".bai",
                     names["BamFilteredPath"], names["BamFilteredPath"] + ".bai", names["MappingQC"]):
            os.replace(work / name, output_dir / name)
    return outputs
