"""
Generate a small synthetic nanopore-like dataset for testing the workflow:
  - ref/transcripts.fa.gz : 40 transcripts from 15 genes (GENCODE-style headers)
  - ref/annotation.gtf.gz : matching 'gene' and 'transcript' features
  - data/*.fastq.gz       : 3 control + 3 treated samples, ~1,000 reads each, ~5% errors;
                            half of the reads are reverse-complemented, as in unstranded PCR-cDNA
Genes G01-G03 are 4x up and G04-G05 4x down in 'treated', so DESeq2 has something to find.

Usage: python make_test_data.py   (run inside the .test/ folder)
"""

import gzip
import random

random.seed(7)
BASES = "ACGT"
N_GENES, READS_PER_SAMPLE, ERROR_RATE = 15, 1000, 0.05


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def rand_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


def mutate(seq):
    out = []
    for b in seq:
        r = random.random()
        if r < ERROR_RATE * 0.4:  # substitution
            out.append(random.choice(BASES.replace(b, "")))
        elif r < ERROR_RATE * 0.7:  # deletion
            continue
        elif r < ERROR_RATE:  # insertion
            out.extend([b, random.choice(BASES)])
        else:
            out.append(b)
    return "".join(out)


# Genes with 2-4 isoforms built from shared exons (so some reads are ambiguous)
transcripts, gtf_lines = [], []
pos = 1000
for g in range(1, N_GENES + 1):
    gene_id, gene_name = f"ENSGTEST{g:05d}.1", f"G{g:02d}"
    exons = [rand_seq(random.randint(250, 600)) for _ in range(5)]
    n_iso = random.randint(2, 4)
    g_start = pos
    for i in range(1, n_iso + 1):
        keep = sorted(random.sample(range(5), random.randint(3, 5)))
        seq = "".join(exons[k] for k in keep)
        tx_id = f"ENSTTEST{g:05d}{i}.1"
        header = f"{tx_id}|{gene_id}|-|-|{gene_name}-20{i}|{gene_name}|{len(seq)}|protein_coding|"
        transcripts.append((tx_id, gene_name, header, seq))
        attrs = (f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_type "protein_coding"; '
                 f'gene_name "{gene_name}";')
        gtf_lines.append(f"chrT\tTEST\ttranscript\t{pos}\t{pos + len(seq)}\t.\t+\t.\t{attrs}")
    gtf_lines.insert(len(gtf_lines) - n_iso,
                     f'chrT\tTEST\tgene\t{g_start}\t{pos + 3000}\t.\t+\t.\tgene_id "{gene_id}"; '
                     f'gene_type "protein_coding"; gene_name "{gene_name}";')
    pos += 5000

with gzip.open("ref/transcripts.fa.gz", "wt") as f:
    for _, _, header, seq in transcripts:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")

with gzip.open("ref/annotation.gtf.gz", "wt") as f:
    f.write("##description: synthetic test annotation\n")
    f.write("\n".join(gtf_lines) + "\n")

base_weight = {tx_id: random.uniform(1, 10) for tx_id, *_ in transcripts}
fold = {"G01": 4, "G02": 4, "G03": 4, "G04": 0.25, "G05": 0.25}

for condition in ("control", "treated"):
    for rep in (1, 2, 3):
        name = f"{condition}_rep{rep}"
        weights = []
        for tx_id, gene_name, _, _ in transcripts:
            w = base_weight[tx_id] * random.uniform(0.8, 1.25)  # biological noise
            if condition == "treated":
                w *= fold.get(gene_name, 1)
            weights.append(w)
        with gzip.open(f"data/{name}.fastq.gz", "wt") as f:
            for r in range(READS_PER_SAMPLE):
                tx_id, _, _, seq = random.choices(transcripts, weights=weights)[0]
                # 3'-biased fragment, typical of cDNA reads that start from the polyA tail
                start = random.randint(0, max(0, len(seq) - 400)) if random.random() < 0.6 else 0
                read = mutate(seq[start:])
                if random.random() < 0.5:  # cDNA reads come from either strand
                    read = revcomp(read)
                qual = "".join(chr(33 + random.randint(8, 20)) for _ in read)
                f.write(f"@{name}_read{r}\n{read}\n+\n{qual}\n")

with open("config/samples.tsv", "w") as f:
    f.write("sample\tcondition\tbatch\tfastq\n")
    for condition in ("control", "treated"):
        for rep in (1, 2, 3):
            batch = "A" if rep == 1 else "B"  # balanced technical batch, as in the SG-NEx design
            f.write(f"{condition}_rep{rep}\t{condition}\t{batch}\tdata/{condition}_rep{rep}.fastq.gz\n")

print(f"{len(transcripts)} transcripts, {N_GENES} genes, 6 samples x {READS_PER_SAMPLE} reads")
