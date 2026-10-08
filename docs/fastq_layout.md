# Inspect original FASTQ layout

Before choosing an extraction pattern, inspect the original, untrimmed inputs:

```bash
python scripts/inspect_fastq_layout.py \
    --manifest samples.tsv \
    --max-reads 100000 \
    --output fastq_layout.json
```

The script uses only the Python standard library, one CPU, and a bounded sample
of the first 100,000 records per unique FASTQ. It supports gzip and plain
four-line FASTQ files. It never trims reads or modifies input files and refuses
to replace an existing report. Relative input paths resolve against the manifest
directory. Shared resolved input paths are reported and inspected once.

The JSON report contains read lengths, frequent initial sequences, per-position
base counts and entropy, and exact barcode-position counts for every manifest
barcode in both forward and reverse-complement orientations. Offsets are
zero-based: offset 8 means that eight bases precede the barcode. The default
search window is the first 80 bases; change it with `--scan-bases` when needed.
At `--barcode-offset` (default 8), it also counts forward exact matches and
matches allowing at most one mismatch. The denominator for displayed percentages
is all examined records; the report also records how many reads are long enough
for the tested barcode position.

This is a layout diagnostic, not a full-library read count or a random sample.
Short motifs can match by chance. A strong barcode peak suggests its position
but does not establish that preceding bases are UMIs. Initial sequence diversity
also does not prove UMI identity: consult the library preparation protocol and
sequencing setup. A barcode used for demultiplexing may be in an index read and
absent from R1. Confirm whether adapter sequences remain after a putative UMI
and barcode before constructing the final extraction pattern.

For already demultiplexed biological samples, correct accidental reuse of a
FASTQ path in the manifest before downstream analysis. Preserve the original
manifest and write the corrected manifest to a new file.
