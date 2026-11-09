"""
Prepare the transcriptome reference:
  1. Write a FASTA whose headers are plain transcript IDs. GENCODE headers look like
     ENST00000456328.2|ENSG00000290825.1|...|DDX11L2-202|...; only the first field is kept.
  2. Build a transcript -> gene table (transcript_id, gene_id, gene_name, gene_type) from the GTF.
"""

import gzip
import re
import sys

sys.stderr = open(snakemake.log[0], "w")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


# 1. FASTA with clean transcript IDs
n_tx = 0
with open_text(snakemake.input.fasta) as fin, open(snakemake.output.fasta, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            tx_id = line[1:].strip().split()[0].split("|")[0]
            fout.write(f">{tx_id}\n")
            n_tx += 1
        else:
            fout.write(line)
print(f"Transcripts written: {n_tx}", file=sys.stderr)

# 2. transcript -> gene table from GTF "transcript" features
attr_re = re.compile(r'(\S+) "([^"]*)"')
rows = {}
with open_text(snakemake.input.gtf) as fin:
    for line in fin:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "transcript":
            continue
        attrs = dict(attr_re.findall(fields[8]))
        tx_id = attrs.get("transcript_id")
        if tx_id is None:
            continue
        rows[tx_id] = (
            attrs.get("gene_id", "NA"),
            attrs.get("gene_name", attrs.get("gene_id", "NA")),
            attrs.get("gene_type", attrs.get("gene_biotype", "NA")),
        )

with open(snakemake.output.tx2gene, "w") as fout:
    fout.write("transcript_id\tgene_id\tgene_name\tgene_type\n")
    for tx_id, (gene_id, gene_name, gene_type) in rows.items():
        fout.write(f"{tx_id}\t{gene_id}\t{gene_name}\t{gene_type}\n")
print(f"Transcript-to-gene pairs: {len(rows)}", file=sys.stderr)

if not rows:
    sys.exit("No 'transcript' features found in the GTF.")
