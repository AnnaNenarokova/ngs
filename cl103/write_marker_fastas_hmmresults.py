#!/usr/bin/env python3
"""
Make one FASTA per TIGRFAM marker (query), containing up to 1 best hit per proteome.

Inputs:
  - A directory of HMMER --tblout files (one per proteome), produced by:
      hmmsearch --noali --cut_ga --tblout PROTEOME.tbl TIGR_MARKERS_ALL.hmm PROTEOME.fasta
  - A directory of proteome FASTA files (same filenames used to name the .tbl files)

Logic:
  1) For each (proteome, marker): choose best hit by highest bit score (tie-break: lowest evalue)
  2) Within each proteome: if the same protein is selected for multiple markers, keep the highest score assignment
  3) Extract sequences from proteome FASTA and write: out_dir/<marker>.fasta
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO


@dataclass(frozen=True)
class Hit:
    proteome: str
    marker: str
    target_id: str
    evalue: float
    score: float


def sanitize_filename(name: str) -> str:
    # keep it filesystem-safe
    name = name.strip()
    name = re.sub(r"[^\w.\-]+", "_", name)
    return name or "marker"


def parse_tblout(tbl_path: Path) -> Iterable[Tuple[str, str, float, float]]:
    """
    Parse HMMER3 --tblout.
    Yields tuples: (target_name, query_name, evalue, score)
    """
    with tbl_path.open("r", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 6:
                continue

            # HMMER tblout columns (typical):
            # target_name target_accession query_name query_accession full_E-value full_score full_bias ...
            target = parts[0]
            query = parts[2]
            try:
                evalue = float(parts[4])
                score = float(parts[5])
            except ValueError:
                continue

            yield target, query, evalue, score


def choose_best_per_marker_per_proteome(tbl_path: Path, proteome_name: str) -> Dict[Tuple[str, str], Hit]:
    """
    From one tblout (one proteome), choose the best hit per marker by highest score.
    Returns dict keyed by (proteome, marker).
    """
    best: Dict[Tuple[str, str], Hit] = {}

    for target, marker, evalue, score in parse_tblout(tbl_path):
        key = (proteome_name, marker)
        h = Hit(proteome=proteome_name, marker=marker, target_id=target, evalue=evalue, score=score)

        if key not in best:
            best[key] = h
            continue

        prev = best[key]
        # best = highest score; tie-break lowest evalue
        if (h.score > prev.score) or (h.score == prev.score and h.evalue < prev.evalue):
            best[key] = h

    return best


def resolve_within_proteome_conflicts(best_hits: Dict[Tuple[str, str], Hit]) -> Dict[Tuple[str, str], Hit]:
    """
    If within the same proteome the same target_id is selected for multiple markers,
    keep the assignment with highest score (tie-break lowest evalue), drop the others.
    """
    by_prot_and_target: Dict[Tuple[str, str], List[Hit]] = defaultdict(list)
    for (_, _), hit in best_hits.items():
        by_prot_and_target[(hit.proteome, hit.target_id)].append(hit)

    kept = dict(best_hits)

    for (proteome, target_id), hits in by_prot_and_target.items():
        if len(hits) <= 1:
            continue
        hits.sort(key=lambda h: (-h.score, h.evalue))
        winner = hits[0]
        for loser in hits[1:]:
            kept.pop((loser.proteome, loser.marker), None)

    return kept


def build_id_index_keys(record_id: str) -> List[str]:
    """
    Try a few reasonable keys for matching HMMER target_name to FASTA record.
    Adjust here if your headers are funky.
    """
    keys = [record_id]
    # If record.id has pipes/spaces, record.id from Biopython is usually first token already,
    # but we keep a couple of fallbacks.
    keys.append(record_id.split()[0])
    keys.append(record_id.split("|")[-1])
    return list(dict.fromkeys(keys))  # preserve order, unique


def extract_sequences(
    proteomes_dir: Path,
    proteome_ext: str,
    selected_hits: Dict[Tuple[str, str], Hit],
) -> Dict[str, List]:
    """
    Returns marker -> list of SeqRecords.
    Loads each proteome once, extracts only needed IDs.
    """
    # proteome -> marker -> target_id
    wanted: Dict[str, Dict[str, str]] = defaultdict(dict)
    for (prot, marker), hit in selected_hits.items():
        wanted[prot][marker] = hit.target_id

    marker_records: Dict[str, List] = defaultdict(list)

    for proteome_name, marker_to_target in wanted.items():
        proteome_path = proteomes_dir / proteome_name
        if not proteome_path.exists() and proteome_ext and not proteome_name.endswith(proteome_ext):
            # try adding extension if user passed bare names
            proteome_path = proteomes_dir / f"{proteome_name}{proteome_ext}"

        if not proteome_path.exists():
            raise FileNotFoundError(f"Missing proteome FASTA for {proteome_name}: {proteome_path}")

        # For fast membership checks
        target_to_marker = {t: m for m, t in marker_to_target.items()}
        needed_targets = set(target_to_marker.keys())

        for rec in SeqIO.parse(str(proteome_path), "fasta"):
            # match using a few common key variants
            matched_target: Optional[str] = None
            for k in build_id_index_keys(rec.id):
                if k in needed_targets:
                    matched_target = k
                    break
            if matched_target is None:
                continue

            marker = target_to_marker[matched_target]

            # Helpful standardized header: proteome|target_id
            rec.id = f"{proteome_name}|{matched_target}"
            rec.description = ""

            marker_records[marker].append(rec)

    return marker_records


def main() -> int:
    ap = argparse.ArgumentParser(
        description="From hmmsearch --tblout (TIGRFAM markers), write one FASTA per marker with up to 1 best hit per proteome."
    )
    ap.add_argument("--tbl_dir", required=True, type=Path, help="Directory of per-proteome hmmsearch --tblout files")
    ap.add_argument("--proteomes_dir", required=True, type=Path, help="Directory of proteome FASTA files")
    ap.add_argument("--out_dir", required=True, type=Path, help="Output directory for marker FASTAs")
    ap.add_argument("--proteome_ext", default=".fasta", help="Proteome FASTA extension (default: .fasta)")
    ap.add_argument("--tbl_ext", default=".tbl", help="tblout file extension (default: .tbl)")
    ap.add_argument("--resolve_conflicts", action="store_true", help="Enforce one marker per protein within each proteome")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Collect best hit per (proteome, marker)
    best_hits: Dict[Tuple[str, str], Hit] = {}

    tbl_files = sorted([p for p in args.tbl_dir.iterdir() if p.is_file() and p.name.endswith(args.tbl_ext)])
    if not tbl_files:
        raise SystemExit(f"No *{args.tbl_ext} files found in {args.tbl_dir}")

    for tbl_path in tbl_files:
        # proteome name is the tbl filename without the .tbl extension
        proteome_name = tbl_path.name[: -len(args.tbl_ext)]
        per_tbl_best = choose_best_per_marker_per_proteome(tbl_path, proteome_name)
        best_hits.update(per_tbl_best)

    # 2) Optional: resolve within-proteome conflicts (same target picked for multiple markers)
    if args.resolve_conflicts:
        best_hits = resolve_within_proteome_conflicts(best_hits)

    # 3) Extract sequences from proteomes and write per marker FASTA
    marker_records = extract_sequences(args.proteomes_dir, args.proteome_ext, best_hits)

    n_written = 0
    for marker, records in marker_records.items():
        out_path = args.out_dir / f"{sanitize_filename(marker)}.fasta"
        SeqIO.write(records, str(out_path), "fasta")
        n_written += 1

    print(f"Selected assignments: {len(best_hits)}")
    print(f"Wrote marker FASTAs:   {n_written} -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
