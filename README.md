# nanopore-rnaseq-snakemake

[![Tests](https://github.com/dlagunasantof/nanopore-rnaseq-snakemake/actions/workflows/ci.yml/badge.svg)](https://github.com/dlagunasantof/nanopore-rnaseq-snakemake/actions/workflows/ci.yml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A reproducible Snakemake workflow for **Oxford Nanopore long-read RNA-seq**: read QC, alignment to the transcriptome, isoform-level quantification and gene-level differential expression. Every step runs in its own pinned Conda environment, and the workflow is tested automatically on a small synthetic dataset with GitHub Actions.

## Workflow

```
FASTQ ──► seqkit (optional subsampling) ──► NanoPlot (read QC) ───────────────┐
  │                                                                           ▼
  └──► minimap2 -ax map-ont (transcriptome) ──► samtools flagstat ──────► MultiQC report
                     │
                     ▼
               NanoCount (EM isoform quantification)
                     │
                     ▼
         transcript / gene count matrices ──► DESeq2 (results, PCA, MA plot)
```

| Step | Tool | Output |
| --- | --- | --- |
| Optional subsampling | seqkit | `results/reads/` |
| Read QC | NanoPlot | `results/qc/nanoplot/` |
| Transcriptome alignment | minimap2 (`-ax map-ont -p 0 -N 10`), samtools | `results/alignments/` |
| Transcript quantification | NanoCount | `results/nanocount/` |
| Count matrices | pandas | `results/counts/` |
| Differential expression | DESeq2 | `results/deseq2/` |
| Summary report | MultiQC | `results/multiqc/multiqc_report.html` |

**Design choices**

- Reads are aligned to the **transcriptome**, keeping secondary alignments, so that NanoCount can resolve reads compatible with several isoforms using an expectation–maximisation algorithm.
- Transcript counts are summed per gene (using the GTF) for DESeq2; transcript-level counts and TPM are kept for isoform analyses.
- DESeq2 runs only when the sample sheet has at least two conditions with two or more replicates each.

## Quick start

### 1. Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda "snakemake>=8"
conda activate snakemake
```

### 2. Test on the bundled synthetic data (about 5 minutes the first time, mostly creating environments)

```bash
git clone https://github.com/dlagunasantof/nanopore-rnaseq-snakemake.git
cd nanopore-rnaseq-snakemake
snakemake -s workflow/Snakefile --directory .test --sdm conda --cores 2
```

The test data (`.test/`) contain 6 samples × 1,000 simulated nanopore-like reads from 40 transcripts; three genes are up-regulated and two down-regulated in `treated`. Regenerate them with `python .test/make_test_data.py` from inside `.test/`.

### 3. Run on real data

**Reference** (GENCODE human transcriptome + annotation):

```bash
bash scripts/download_reference.sh 47   # check the latest release at gencodegenes.org
```

**Reads**: the example design uses public data from the [Singapore Nanopore Expression Project (SG-NEx)](https://registry.opendata.aws/sgnex/), PCR-cDNA runs of two cell lines (HepG2, K562) in three replicates each:

```bash
conda install -c conda-forge awscli
bash scripts/download_sgnex.sh list HepG2 cDNA          # find the FASTQ files
bash scripts/download_sgnex.sh get <s3-path> data/HepG2_rep1.fastq.gz
```

Then edit `config/samples.tsv` (columns: `sample`, `condition`, `fastq`) and `config/config.yaml`, and run:

```bash
snakemake -s workflow/Snakefile --sdm conda --cores 8
```

By default each sample is subsampled to 500,000 reads (`subsample_reads` in the config) so the full run fits on a modest server; set it to `0` to use all reads.

### Running on a cluster

An example SLURM profile is provided in `profiles/slurm/` (requires `snakemake-executor-plugin-slurm`; set your partition and account first):

```bash
snakemake -s workflow/Snakefile --profile profiles/slurm
```

## Repository structure

```
├── config/                 config.yaml, samples.tsv
├── workflow/
│   ├── Snakefile
│   ├── rules/              one .smk file per stage
│   ├── envs/               pinned Conda environments
│   └── scripts/            reference preparation, count merging, DESeq2
├── profiles/slurm/         example cluster profile
├── scripts/                data download helpers
├── .test/                  synthetic test dataset and config
└── .github/workflows/      continuous integration
```

## Results on SG-NEx data

<!-- After running on real data, add here: number of reads per sample, mapping rate,
     the PCA plot and the number of differentially expressed genes. -->

## Data and tools

- SG-NEx: Chen Y. *et al.* A systematic benchmark of Nanopore long-read RNA sequencing for transcript-level analysis in human cell lines. *Nature Methods* 22, 801–812 (2025).
- [minimap2](https://github.com/lh3/minimap2) · [samtools](https://github.com/samtools/samtools) · [NanoCount](https://github.com/a-slide/NanoCount) · [NanoPlot](https://github.com/wdecoster/NanoPlot) · [seqkit](https://github.com/shenwei356/seqkit) · [DESeq2](https://bioconductor.org/packages/DESeq2) · [MultiQC](https://multiqc.info) · [Snakemake](https://snakemake.github.io)

## Author

Daniela Laguna Santofimio · [GitHub](https://github.com/dlagunasantof) · [LinkedIn](https://linkedin.com/in/daniela-laguna-santofimio)
