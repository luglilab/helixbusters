"""Real BAM fixtures, independent of external mapping programs."""
import pysam

HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"},
          "SQ": [{"SN": name, "LN": 1000}
                 for name in ("chr1", "chrX", "chrY", "chrM", "scaffold_1")]}


def read(start=100, umi="AAAAAAAA", flag=0, cigar="20M", ref=0,
         name=None, mapq=60, rx=None):
    record = pysam.AlignedSegment()
    record.query_name = name if name is not None else f"instrument_lane_read_{umi}"
    record.flag = flag
    record.reference_id = ref
    record.reference_start = start
    record.mapping_quality = mapq
    record.cigarstring = cigar
    length = sum(n for op, n in record.cigartuples or [] if op in {0, 1, 4, 7, 8})
    record.query_sequence = "A" * (length or 20)
    record.query_qualities = pysam.qualitystring_to_array("I" * (length or 20))
    if rx is not None:
        record.set_tag("RX", rx)
    return record


def write_bam(path, records, sort=True, header=None):
    if sort:
        records = sorted(records, key=lambda r: (
            r.reference_id if r.reference_id >= 0 else 10**9, r.reference_start))
    with pysam.AlignmentFile(str(path), "wb", header=header or HEADER) as out:
        for record in records:
            out.write(record)
    return path
