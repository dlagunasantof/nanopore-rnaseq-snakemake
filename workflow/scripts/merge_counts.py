"""
Merge oarfish outputs (one per sample) into:
  - transcript_counts.tsv : estimated read counts (num_reads), transcripts x samples
  - transcript_cpm.tsv    : counts per million reads, transcripts x samples
                            (no length correction: a long read represents one whole molecule)
  - gene_counts.tsv       : estimated counts summed per gene, rounded for DESeq2
"""

import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

tx2gene = pd.read_csv(snakemake.input.tx2gene, sep="\t", dtype=str)
# Match on version-less IDs so ENST0000123.4 in the FASTA matches ENST0000123.4 or ENST0000123 in the GTF
tx2gene["tx_key"] = tx2gene["transcript_id"].str.replace(r"\.\d+$", "", regex=True)
tx2gene = tx2gene.drop_duplicates("tx_key").set_index("tx_key")

counts = {}
for sample, path in zip(snakemake.params.samples, snakemake.input.counts):
    df = pd.read_csv(path, sep="\t").set_index("tname")
    counts[sample] = df["num_reads"]

tx_counts = pd.DataFrame(counts).fillna(0.0)
tx_counts = tx_counts[tx_counts.sum(axis=1) > 0]  # drop transcripts with no reads in any sample
tx_cpm = tx_counts / tx_counts.sum(axis=0) * 1e6
tx_counts.index.name = tx_cpm.index.name = "transcript_id"

tx_counts.round(3).to_csv(snakemake.output.tx_counts, sep="\t")
tx_cpm.round(3).to_csv(snakemake.output.tx_cpm, sep="\t")

# Gene level
keys = tx_counts.index.str.replace(r"\.\d+$", "", regex=True)
mapped = tx2gene.reindex(keys)
n_unmapped = mapped["gene_id"].isna().sum()
if n_unmapped:
    print(f"Warning: {n_unmapped} transcripts have no gene in the GTF and are skipped.", file=sys.stderr)

gene_counts = tx_counts.copy()
gene_counts["gene_id"] = mapped["gene_id"].values
gene_counts["gene_name"] = mapped["gene_name"].values
gene_counts = gene_counts.dropna(subset=["gene_id"])
gene_counts = gene_counts.groupby(["gene_id", "gene_name"]).sum().round().astype(int)
gene_counts.to_csv(snakemake.output.gene_counts, sep="\t")

print(
    f"Transcripts: {tx_counts.shape[0]}, genes: {gene_counts.shape[0]}, samples: {tx_counts.shape[1]}",
    file=sys.stderr,
)
