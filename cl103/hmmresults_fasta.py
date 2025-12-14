import os
from pathlib import Path
from Bio import SeqIO

def prepare_fastas(
    hmm_report_dir,
    proteomes_dir,
    cog_dir,
    result_dir,
    hmm_ext=".hmm",
    proteome_ext=".fasta",
    n_best=1,
    max_evalue=1e-7,
):
    """
    For each COG fasta in `cog_dir`, append matching HMM-hit sequences
    from proteomes in `proteomes_dir`, and write the combined fasta to `result_dir`.

    Expected structure of `prepare_hmm_dict` output:
      hmm_dict[proteome][cog] is a dict/set-like container keyed by sseqid.
      This function replaces values with SeqRecord objects when found.
    """
    hmm_report_dir = Path(hmm_report_dir)
    proteomes_dir = Path(proteomes_dir)
    cog_dir = Path(cog_dir)
    result_dir = Path(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)

    print("Parsing hmm_reports")
    hmm_dict = prepare_hmm_dict(
        str(hmm_report_dir),
        n_best=n_best,
        max_evalue=max_evalue,
        hmm_ext=hmm_ext,
    )

    print("Parsing proteome sequences")
    for proteome_file in listdir_nohidden(str(proteomes_dir)):
        if proteome_ext and not proteome_file.endswith(proteome_ext):
            continue

        proteome_path = proteomes_dir / proteome_file
        proteome_hits = hmm_dict.get(proteome_file)
        if not proteome_hits:
            # No HMM report (or empty) for this proteome; skip safely.
            continue

        # Precompute all IDs we need from this proteome to avoid looping over all COGs per record
        needed_ids = set()
        for cog, hits in proteome_hits.items():
            needed_ids.update(hits.keys() if hasattr(hits, "keys") else hits)

        for record in SeqIO.parse(str(proteome_path), "fasta"):
            sseqid = record.id
            if sseqid not in needed_ids:
                continue

            # Populate record for all COGs that expect it
            for cog, hits in proteome_hits.items():
                if sseqid in hits:
                    hits[sseqid] = record

    print("Writing down sequences")
    for cog_file in listdir_nohidden(str(cog_dir)):
        cog_path = cog_dir / cog_file
        out_path = result_dir / cog_file

        out_records = list(SeqIO.parse(str(cog_path), "fasta"))

        for proteome_file, proteome_hits in hmm_dict.items():
            proteome_cog_dict = proteome_hits.get(cog_file)
            if not proteome_cog_dict:
                # Optional: comment out if too noisy
                # print(f"{cog_file} not in hmm_
