#!/usr/bin/env python3
"""Fetch a UCSC iGenome BWA index and write a Helixbusters reference catalog."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import sys
import tarfile
import subprocess
import time


BASE_URL = "https://s3.amazonaws.com/igenomes.illumina.com"
BUILD_INFO = {
    "hg19": ("Homo_sapiens", "UCSC", "hg19", "Homo_sapiens_UCSC_hg19.tar.gz"),
    "hg38": ("Homo_sapiens", "UCSC", "hg38", "Homo_sapiens_UCSC_hg38.tar.gz"),
    "mm10": ("Mus_musculus", "UCSC", "mm10", "Mus_musculus_UCSC_mm10.tar.gz"),
}
INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
BOYLE_BLACKLISTS = {
    "hg19": "hg19-blacklist.v2.bed.gz",
    "hg38": "hg38-blacklist.v2.bed.gz",
    "mm10": "mm10-blacklist.v2.bed.gz",
}
BLACKLIST_BASE_URL = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists"


def download(url, destination, rate_limit=None, attempts=8):
    """Resume a persistent partial file with curl; retain it after failure."""
    if not shutil.which("curl"):
        raise ValueError("curl is required for resumable downloads; load it in your environment")
    command = ["curl", "--fail", "--location", "--show-error", "--silent",
               "--connect-timeout", "30", "--speed-limit", "1", "--speed-time", "120",
               "--continue-at", "-", "--output", str(destination)]
    if rate_limit:
        command.extend(["--limit-rate", rate_limit])
    command.append(url)
    for attempt in range(1, attempts + 1):
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode == 0:
            if not destination.is_file() or destination.stat().st_size == 0:
                raise ValueError(f"Empty download: {url}")
            digest = hashlib.sha256()
            with destination.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            return digest.hexdigest()
        size = destination.stat().st_size if destination.exists() else 0
        detail = result.stderr.strip()
        # Invalid options, certificates, HTTP errors and unsupported ranges
        # require operator intervention rather than repeated identical requests.
        if result.returncode in (2, 3, 22, 33, 60) or attempt == attempts:
            raise ValueError(f"Download failed at {size} bytes: {detail}; partial retained at {destination}")
        delay = min(2 ** attempt, 60)
        print(f"Download interrupted at {size} bytes; retry {attempt + 1}/{attempts} "
              f"in {delay}s: {detail}", flush=True)
        time.sleep(delay)


def find_bwa_prefix(directory):
    """Find a nonempty classic BWA index without modifying installed references."""
    candidates = [directory / "genome.fa"]
    candidates.extend(Path(str(path)[:-4]) for path in sorted(directory.glob("*/genome.fa.amb")))
    for prefix in candidates:
        if all(Path(str(prefix) + suffix).is_file() and
               Path(str(prefix) + suffix).stat().st_size > 0 for suffix in INDEX_SUFFIXES):
            return prefix
    raise ValueError(f"No complete nonempty classic BWA index under {directory}")


def extract_bwa_index(archive, genome_root, target_dir):
    """Extract only the five classic BWA index files from the iGenome tarball."""
    relative_root = f"{genome_root}/Sequence/BWAIndex/"
    with tarfile.open(archive, "r:gz") as tar:
        members = [member for member in tar.getmembers()
                   if member.isfile() and member.name.lstrip("./").startswith(relative_root)]
        prefixes = set()
        for member in members:
            name = member.name.lstrip("./")
            if name.endswith(".amb"):
                prefixes.add(name[:-4])
        complete = [prefix for prefix in prefixes
                    if all(any(m.name.lstrip("./") == prefix + suffix for m in members)
                           for suffix in INDEX_SUFFIXES)]
        if not complete:
            raise ValueError(f"No complete classic BWA index found in iGenome archive under {relative_root}")
        prefix = sorted(complete)[-1]
        target_dir.mkdir(parents=True, exist_ok=True)
        for suffix in INDEX_SUFFIXES:
            member = next(m for m in members if m.name.lstrip("./") == prefix + suffix)
            source = tar.extractfile(member)
            if source is None:
                raise ValueError(f"Could not read {member.name} from iGenome archive")
            with source, (target_dir / ("genome.fa" + suffix)).open("wb") as output:
                shutil.copyfileobj(source, output)
        version = PurePosixPath(prefix).parent.name
    return version


def validate_blacklist(path):
    """Read and validate BED columns and gzip integrity before using a list."""
    intervals = 0
    with path.open("rb") as raw:
        magic = raw.read(2)
    named_gzip = path.name.lower().endswith(".gz")
    if named_gzip and magic != b"\x1f\x8b":
        raise ValueError(f"{path} has a .gz suffix but is not gzip-compressed")
    opener = gzip.open if magic == b"\x1f\x8b" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"Invalid blacklist BED row {number}: fewer than 3 columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"Invalid blacklist BED coordinates on row {number}") from error
            if start < 0 or end <= start:
                raise ValueError(f"Invalid blacklist BED interval on row {number}")
            intervals += 1
    if not intervals:
        raise ValueError("Blacklist BED contains no intervals")


def prepare_blacklist(genome, cache, rate_limit=None, attempts=8):
    filename = BOYLE_BLACKLISTS.get(genome)
    if filename is None:
        raise ValueError(f"Boyle-Lab has no listed blacklist for {genome}; supply a build-matched BED manually")
    directory = cache / "blacklists"
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / filename
    url = f"{BLACKLIST_BASE_URL}/{filename}"
    if not destination.is_file() or destination.stat().st_size == 0:
        partial = destination.with_suffix(destination.suffix + ".partial")
        print(f"Downloading {url}", flush=True)
        digest = download(url, partial, rate_limit, attempts)
        validate_blacklist(partial)
        partial.replace(destination)
    else:
        validate_blacklist(destination)
        digest = hashlib.sha256(destination.read_bytes()).hexdigest()
    provenance = {"source": "Boyle-Lab/Blacklist", "url": url, "genome": genome,
                  "file": str(destination), "sha256": digest}
    (directory / (filename + ".provenance.json")).write_text(json.dumps(provenance, indent=2) + "\n")
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True, choices=sorted(BUILD_INFO),
                        help="Supported iGenomes UCSC build: hg19, hg38 or mm10 (not mm39)")
    parser.add_argument("--blacklist", help="Optional existing build-specific BED; otherwise download the Boyle-Lab list")
    parser.add_argument("--blacklist-genome",
                        help="Required with --blacklist; explicit build label, which must match --genome")
    parser.add_argument("--cache-dir", required=True, help="Persistent reference cache directory on the HPC")
    parser.add_argument("--output-config", required=True, help="Output references.json path")
    parser.add_argument("--rate-limit", help="curl transfer limit, e.g. 100K or 10M (bytes/s)")
    parser.add_argument("--attempts", type=int, default=8, help="Maximum download attempts (default: 8)")
    parser.add_argument("--igenomes-base", help="Existing local iGenomes root; do not download the archive")
    parser.add_argument("--base-url", default=BASE_URL, help="iGenomes archive mirror base URL")
    args = parser.parse_args()
    if args.attempts < 1:
        parser.error("--attempts must be positive")
    config_path = Path(args.output_config).expanduser().resolve()
    if config_path.exists():
        parser.error(f"output catalog already exists: {config_path}; choose a new --output-config")

    genome = args.genome
    cache = Path(args.cache_dir).expanduser().resolve()
    cache.mkdir(parents=True, exist_ok=True)
    if args.blacklist:
        if not args.blacklist_genome:
            parser.error("--blacklist-genome is required when --blacklist is supplied")
        if args.blacklist_genome.lower() != genome:
            parser.error("--blacklist-genome must exactly match --genome for this UCSC iGenome mode")
        blacklist = Path(args.blacklist).expanduser().resolve()
        if not blacklist.is_file() or blacklist.stat().st_size == 0:
            parser.error(f"blacklist does not exist or is empty: {blacklist}")
        try:
            validate_blacklist(blacklist)
        except (OSError, EOFError, ValueError) as error:
            parser.error(f"invalid blacklist: {error}")
    else:
        try:
            blacklist = prepare_blacklist(genome, cache, args.rate_limit, args.attempts)
        except (OSError, EOFError, ValueError) as error:
            parser.error(str(error))

    organism, source, build, archive_name = BUILD_INFO[genome]
    genome_root = f"{organism}/{source}/{build}"
    url = f"{args.base_url.rstrip('/')}/{genome_root}/{archive_name}"
    index_dir = cache / genome / "BWAIndex"
    index_prefix = index_dir / "genome.fa"
    config_path = Path(args.output_config).expanduser().resolve()
    cache.mkdir(parents=True, exist_ok=True)

    provenance_path = cache / genome / "igenome.provenance.json"
    if args.igenomes_base:
        directory = Path(args.igenomes_base).expanduser().resolve() / genome_root / "Sequence/BWAIndex"
        index_prefix = find_bwa_prefix(directory)
        provenance = {"source": "Existing iGenomes installation", "build": genome,
                      "index": str(index_prefix), "archive_sha256": None}
    elif all(Path(str(index_prefix) + suffix).is_file() and
             Path(str(index_prefix) + suffix).stat().st_size > 0 for suffix in INDEX_SUFFIXES):
        provenance = json.loads(provenance_path.read_text()) if provenance_path.is_file() else {
            "source": "Illumina iGenomes", "build": genome, "index": str(index_prefix),
            "archive_sha256": None, "note": "Reused existing BWA index cache; source archive checksum unavailable"}
    else:
        if any(Path(str(index_prefix) + suffix).exists() for suffix in INDEX_SUFFIXES):
            raise ValueError(f"Incomplete existing index at {index_prefix}; choose a new cache directory")
        archive_dir = cache / "downloads"
        archive_dir.mkdir(parents=True, exist_ok=True)
        archive = archive_dir / (archive_name + ".partial")
        print(f"Downloading {url} (partial cache: {archive})", flush=True)
        archive_sha256 = download(url, archive, args.rate_limit, args.attempts)
        version = extract_bwa_index(archive, genome_root, index_dir)
        provenance = {"source": "Illumina iGenomes", "url": url, "build": genome,
                      "index": str(index_prefix), "bwa_index_version": version,
                      "archive_sha256": archive_sha256}
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")

    config = {"schema_version": 1, "references": {genome: {
        "bwa_index": str(index_prefix),
        "blacklist": {"genome": genome, "path": str(blacklist)},
    }}}
    config_path.parent.mkdir(parents=True, exist_ok=True)
    temp_config = config_path.with_suffix(config_path.suffix + ".tmp")
    temp_config.write_text(json.dumps(config, indent=2) + "\n")
    temp_config.replace(config_path)
    print(f"BWA index ready: {index_prefix}")
    print(f"Reference catalog written: {config_path}")
    print(f"Download provenance: {provenance_path}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, EOFError, tarfile.TarError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
