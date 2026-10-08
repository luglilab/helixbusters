"""Build selection, canonical nuclear chromosomes and BED blacklist filtering."""

from bisect import bisect_right
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path


# Canonical nuclear lengths from UCSC goldenPath/<build>/bigZips/<build>.chrom.sizes.
# Order: autosomes, X, Y. Mitochondrial lengths are deliberately not compared:
# hg19 and some GRCh37 distributions use different mitochondrial references.
_LENGTHS = {
    "hg19": [249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
             159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
             115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
             59128983, 63025520, 48129895, 51304566, 155270560, 59373566],
    "hg38": [248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
             159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
             114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
             58617616, 64444167, 46709983, 50818468, 156040895, 57227415],
    "mm10": [195471971, 182113224, 160039680, 156508116, 151834684, 149736546,
             145441459, 129401213, 124595110, 130694993, 122082543, 120129022,
             120421639, 124902244, 104043685, 98207768, 94987271, 90702639,
             61431566, 171031299, 91744698],
    "mm39": [195154279, 181755017, 159745316, 156860686, 151758149, 149588044,
             144995196, 130127694, 124359700, 130530862, 121973369, 120092757,
             120883175, 125139656, 104073951, 98008968, 95294699, 90720763,
             61420004, 169476592, 91455967],
}
_ALIASES = {"hg19": "hg19", "grch37": "hg19", "hg38": "hg38", "grch38": "hg38",
            "mm10": "mm10", "grcm38": "mm10", "mm39": "mm39", "grcm39": "mm39"}
MITOCHONDRIAL = frozenset({"chrM", "chrMT", "MT", "M"})


def normalize_genome(genome):
    if not isinstance(genome, str) or genome.lower() not in _ALIASES:
        raise ValueError("genome must be hg19/GRCh37, hg38/GRCh38, mm10/GRCm38 or mm39/GRCm39")
    return _ALIASES[genome.lower()]


def genome_species(genome):
    return "human" if normalize_genome(genome).startswith("hg") else "mouse"


def canonical_lengths(genome):
    genome = normalize_genome(genome)
    autosomes = 22 if genome_species(genome) == "human" else 19
    names = [f"chr{i}" for i in range(1, autosomes + 1)] + ["chrX", "chrY"]
    return dict(zip(names, _LENGTHS[genome]))


class GenomeFilter:
    """Canonical nuclear-only filter with an explicitly build-labelled BED.

    BED coordinates are zero-based and half-open. An overlap of >=1 aligned
    base excludes the read. D/N gaps and clipped bases are not aligned blocks.
    Input chromosome names are never rewritten in BAM or deduplication outputs.
    """

    def __init__(self, genome, blacklist_bed, blacklist_genome):
        self.genome = normalize_genome(genome)
        if blacklist_genome is None:
            raise ValueError("Declare blacklist_genome explicitly; a BED alone does not identify its build")
        if normalize_genome(blacklist_genome) != self.genome:
            raise ValueError("Blacklist genome does not match the selected genome")
        if blacklist_bed is None:
            raise ValueError("A build-matched blacklist_bed is required when selecting a genome")
        self.lengths = canonical_lengths(self.genome)
        self.aliases = {alias: chrom for chrom in self.lengths for alias in (chrom, chrom[3:])}
        self.path = Path(blacklist_bed).expanduser().resolve()
        raw = self.path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        text = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode("utf-8-sig")
        intervals = defaultdict(list)
        self.input_intervals = self.ignored_intervals = 0
        for number, line in enumerate(text.splitlines(), 1):
            line = line.strip()
            if not line or line.startswith("#") or line.split()[0] in {"track", "browser"}:
                continue
            fields = line.split()
            try:
                start, end = int(fields[1]), int(fields[2])
            except (IndexError, ValueError) as error:
                raise ValueError(f"Invalid BED row at {self.path}:{number}") from error
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval at {self.path}:{number}: require 0 <= start < end")
            self.input_intervals += 1
            chrom = self.aliases.get(fields[0])
            if chrom is None:
                self.ignored_intervals += 1
                continue
            if end > self.lengths[chrom]:
                raise ValueError(f"Blacklist interval exceeds {self.genome} {chrom} at row {number}")
            intervals[chrom].append((start, end))
        if not intervals:
            raise ValueError("Blacklist has no intervals on canonical chromosomes; check file and chromosome names")
        self.intervals = {}
        for chrom, regions in intervals.items():
            merged = []
            for start, end in sorted(regions):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.intervals[chrom] = ([start for start, end in merged], [end for start, end in merged])

    def validate_header(self, references, lengths):
        """Detect a wrong build or naming mismatch before producing filtered BAMs.

        Requires the full canonical nuclear complement and expected lengths.
        Length agreement is a compatibility check, not sequence authentication.
        """
        found = {}
        for name, length in zip(references, lengths):
            chrom = self.aliases.get(name)
            if chrom is None:
                continue
            if chrom in found:
                raise ValueError(f"Ambiguous BAM naming: multiple aliases for {chrom}")
            if length != self.lengths[chrom]:
                raise ValueError(f"BAM/index genome mismatch: {name} length {length} is not {self.genome} ({self.lengths[chrom]})")
            found[chrom] = name
        missing = self.lengths.keys() - found.keys()
        if missing:
            raise ValueError(f"Missing canonical chromosomes for {self.genome}: {', '.join(sorted(missing))}; use UCSC or Ensembl chromosome names")

    def overlaps(self, chrom, start, end):
        regions = self.intervals.get(self.aliases.get(chrom))
        if regions is None:
            return False
        starts, ends = regions
        index = bisect_right(ends, start)
        return index < len(starts) and starts[index] < end

    def exclusion_reason(self, read):
        chrom = read.reference_name
        if chrom in MITOCHONDRIAL:
            return "mitochondrial"
        if chrom not in self.aliases:
            return "noncanonical"
        if any(self.overlaps(chrom, start, end) for start, end in read.get_blocks()):
            return "blacklist"
        return None

    def metadata(self):
        return {"genome": self.genome, "species": genome_species(self.genome),
                "canonical_nuclear_chromosomes": list(self.lengths),
                "exclude_mitochondrial": True, "blacklist_path": str(self.path),
                "blacklist_genome": self.genome, "blacklist_sha256": self.sha256,
                "blacklist_input_intervals": self.input_intervals,
                "blacklist_noncanonical_intervals_ignored": self.ignored_intervals,
                "blacklist_merged_intervals": sum(len(starts) for starts, ends in self.intervals.values()),
                "overlap_policy": "at least one aligned reference base; CIGAR blocks; per read"}


def make_genome_filter(genome=None, blacklist_bed=None, blacklist_genome=None):
    if genome is None:
        if blacklist_bed is not None or blacklist_genome is not None:
            raise ValueError("Select genome before supplying a blacklist")
        return None  # Backward-compatible explicitly unconfigured mode.
    return GenomeFilter(genome, blacklist_bed, blacklist_genome)


def load_reference_config(config_path, genome, aligner="bwa"):
    """Select local index and build-labelled blacklist; resolve relative paths."""
    genome = normalize_genome(genome)
    if aligner not in {"bwa", "bowtie2"}:
        raise ValueError("aligner must be 'bwa' or 'bowtie2'")
    config_path = Path(config_path).expanduser().resolve()
    config = json.loads(config_path.read_text())
    if config.get("schema_version") != 1:
        raise ValueError("Reference config requires schema_version=1")
    profiles = {}
    for name, entry in config.get("references", {}).items():
        normalized = normalize_genome(name)
        if normalized in profiles:
            raise ValueError(f"Duplicate reference aliases for {normalized}")
        profiles[normalized] = entry
    if genome not in profiles:
        raise ValueError(f"Reference config has no entry for {genome}")
    entry = profiles[genome]
    blacklist = entry.get("blacklist", {})
    if not isinstance(blacklist, dict):
        raise ValueError("blacklist must contain path and genome")

    def resolve(value):
        if not isinstance(value, str) or not value:
            raise ValueError(f"Missing reference path for {genome}/{aligner}")
        path = Path(value).expanduser()
        return str((config_path.parent / path).resolve()) if not path.is_absolute() else str(path)

    return {"genome": genome, "genome_index": resolve(entry.get(f"{aligner}_index")),
            "blacklist_bed": resolve(blacklist.get("path")),
            "blacklist_genome": blacklist.get("genome")}
