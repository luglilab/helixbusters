#!/usr/bin/env python3
"""Validate a Helixbusters single-end samplesheet and write a TSV manifest."""

import argparse
from pathlib import Path
import sys

import pandas as pd


REQUIRED = {"Sample", "Replicate", "Group", "PathReadForward", "SampleBarcodeForward"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samplesheet", help="Input .xlsx, .xls or .csv")
    parser.add_argument("manifest", help="Output tab-separated manifest")
    parser.add_argument("--check-fastq", action="store_true",
                        help="Require FASTQ files to exist on this machine")
    args = parser.parse_args()

    sheet_path = Path(args.samplesheet).expanduser().resolve()
    suffix = sheet_path.suffix.lower()
    if suffix == ".csv":
        frame = pd.read_csv(sheet_path)
    elif suffix in {".xls", ".xlsx"}:
        frame = pd.read_excel(sheet_path)
    else:
        parser.error("samplesheet must be .csv, .xls or .xlsx")

    missing = REQUIRED - set(frame.columns)
    if missing:
        parser.error(f"missing required columns: {', '.join(sorted(missing))}")
    paired_columns = {"PathReadReverse", "SampleBarcodeReverse"}
    if paired_columns & set(frame.columns):
        parser.error("this first Nextflow workflow supports single-end samplesheets only")
    if frame.empty:
        parser.error("samplesheet has no sample rows")
    if frame["Sample"].isna().any() or frame["Sample"].duplicated().any():
        parser.error("Sample values must be present and unique")
    if frame[["Replicate", "Group", "PathReadForward", "SampleBarcodeForward"]].isna().any().any():
        parser.error("Replicate, Group, FASTQ path and barcode must be present on every row")

    rows = []
    for record in frame.to_dict(orient="records"):
        sample = str(record["Sample"]).strip()
        if not sample or sample in {".", ".."} or any(c in sample for c in "/\\\t\r\n"):
            parser.error(f"invalid sample name: {sample!r}")
        barcode = str(record["SampleBarcodeForward"]).strip().upper()
        if not barcode or set(barcode) - set("ACGT"):
            parser.error(f"SampleBarcodeForward for {sample} must contain only A/C/G/T")
        fastq = Path(str(record["PathReadForward"])).expanduser()
        if not fastq.is_absolute():
            fastq = (sheet_path.parent / fastq).resolve()
        if args.check_fastq and not fastq.is_file():
            parser.error(f"FASTQ does not exist for {sample}: {fastq}")
        rows.append({
            "sample": sample,
            "replicate": str(record["Replicate"]).strip(),
            "group": str(record["Group"]).strip(),
            "barcode": barcode,
            "fastq": str(fastq.resolve()),
        })

    out = Path(args.manifest).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Validated {len(rows)} single-end samples; manifest: {out}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, ImportError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
