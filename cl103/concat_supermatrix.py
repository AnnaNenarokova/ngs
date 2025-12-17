#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from collections import defaultdict, OrderedDict
from typing import Dict, List, Tuple

from Bio import SeqIO


def get_taxon_id(seq_id: str, mode: str) -> str:
    """
    mode:
      - 'before_pipe'  : taxon = part before first '|'
      - 'full'         : taxon = full seq_id
      - 'first_token'  : taxon = first whitespace-delimited token
    """
    if mode == "before_pipe":
        return seq_id.split("|", 1)[0]
    if mode == "first_token":
        return seq_id.split()[0]
    return seq_id


def read_alignment_fasta(path: Path) -> Dict[str, str]:
    """
    Returns dict: seq_id -> aligned sequence string (uppercase).
    """
    aln = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        aln[rec.id] = str(rec.seq).upper()
    return aln


def ensure_equal_lengths(aln: Dict[str, str], filename: str) -> int:
    lens = {len(s) for s in aln.values()}
    if not lens:
        raise ValueError(f"{filename}: alignment is empty")
    if len(lens) != 1:
        raise ValueError(f"{filename}: sequences have different lengths: {sorted(lens)}")
    return next(iter(lens))


def main() -> int:
    ap = argparse.ArgumentParser(description="Concatenate per-gene MSAs into a supermatrix (FASTA) with optional partitions.")
    ap.add_argument("--msa_dir", required=True, type=Path, help="Directory with per-marker MSAs (FASTA format).")
    ap.add_argument("--out_fasta", required=True, type=Path, help="Output concatenated supermatrix FASTA.")
    ap.add_argument("--msa_ext", default=".aln", help="Only read files ending with this extension (default: .aln).")
    ap.add_argument("--taxon_mode", default="before_pipe", choices=["before_pipe", "full", "first_token"], help="How to extract taxon ID from sequence header.")
    ap.add_argument("--partitions", default=None, type=Path, help="Optional output partitions file (IQ-TREE/RAxML style).")
    ap.add_argument("--sort_taxa", action="store_true", help="Sort taxa alphabetically (otherwise preserve discovery order).")
    ap.add_argument("--sort_markers", action="store_true", help="Sort marker files alphabetically (otherwise filesystem order).")
    args = ap.parse_args()

    msa_files = [p for p in args.msa_dir.iterdir() if p.is_file() and p.name.endswith(args.msa_ext)]
    if not msa_files:
        raise SystemExit(f"No MSA files with extension '{args.msa_ext}' found in: {args.msa_dir}")

    msa_files = sorted(msa_files) if args.sort_markers else list(msa_files)

    # marker -> (taxon -> seq), and marker lengths
    marker_taxon_seq: List[Tuple[str, Dict[str, str], int]] = []
    all_taxa_order = OrderedDict()  # preserves discovery order

    for msa_path in msa_files:
        marker = msa_path.stem  # file name without extension
        aln = read_alignment_fasta(msa_path)
        aln_len = ensure_equal_lengths(aln, msa_path.name)

        taxon_to_seq: Dict[str, str] = {}
        for seq_id, seq in aln.items():
            taxon = get_taxon_id(seq_id, args.taxon_mode)
            if taxon in taxon_to_seq:
                raise ValueError(f"{msa_path.name}: duplicate taxon '{taxon}' after applying taxon_mode={args.taxon_mode}")
            taxon_to_seq[taxon] = seq
            all_taxa_order.setdefault(taxon, None)

        marker_taxon_seq.append((marker, taxon_to_seq, aln_len))

    taxa = list(all_taxa_order.keys())
    if args.sort_taxa:
        taxa = sorted(taxa)

    # Build concatenated sequences
    concat: Dict[str, List[str]] = {t: [] for t in taxa}
    partitions: List[Tuple[str, int, int]] = []

    pos = 1  # 1-based inclusive coordinates for partitions
    for marker, taxon_to_seq, aln_len in marker_taxon_seq:
        start = pos
        end = pos + aln_len - 1
        partitions.append((marker, start, end))
        pos = end + 1

        gap_block = "-" * aln_len
        for taxon in taxa:
            concat[taxon].append(taxon_to_seq.get(taxon, gap_block))

    # Write supermatrix FASTA
    args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.out_fasta.open("w") as out:
        for taxon in taxa:
            out.write(f">{taxon}\n{''.join(concat[taxon])}\n")

    # Write partitions file (optional)
    if args.partitions is not None:
        args.partitions.parent.mkdir(parents=True, exist_ok=True)
        with args.partitions.open("w") as outp:
            for marker, start, end in partitions:
                # IQ-TREE / RAxML-style
                outp.write(f"{marker} = {start}-{end}\n")

    print(f"Wrote supermatrix: {args.out_fasta}")
    if args.partitions is not None:
        print(f"Wrote partitions: {args.partitions}")
    print(f"Taxa: {len(taxa)}; Markers: {len(marker_taxon_seq)}; Total sites: {pos-1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
