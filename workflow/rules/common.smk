import os

import pandas as pd


# ---------------------------------------------------------------------------
# Sample sheet
# ---------------------------------------------------------------------------
samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").set_index(
    "sample", drop=False
)

_required = {"sample", "condition", "fastq"}
_missing = _required - set(samples.columns)
if _missing:
    raise ValueError(f"Sample sheet is missing columns: {', '.join(sorted(_missing))}")
if samples.index.duplicated().any():
    raise ValueError("Sample names in the sample sheet must be unique.")

SAMPLES = samples.index.tolist()

# Technical covariates for DESeq2 (e.g. sequencing kit or batch) must be sample sheet columns
COVARIATES = config.get("deseq2", {}).get("covariates", []) or []
_missing_cov = [c for c in COVARIATES if c not in samples.columns]
if _missing_cov:
    raise ValueError(
        f"deseq2.covariates {_missing_cov} are not columns of the sample sheet {config['samples']}."
    )


wildcard_constraints:
    sample="|".join(SAMPLES),


def get_fastq(wildcards):
    return samples.loc[wildcards.sample, "fastq"]


# ---------------------------------------------------------------------------
# Differential expression is only run when the design allows it:
# at least two conditions, each with at least two replicates.
# ---------------------------------------------------------------------------
def run_deseq2():
    if not config.get("deseq2", {}).get("enabled", True):
        return False
    reps = samples["condition"].value_counts()
    return len(reps) >= 2 and (reps >= 2).all()


def deseq2_outputs():
    if not run_deseq2():
        return []
    return [
        "results/deseq2/deseq2_results.tsv",
        "results/deseq2/normalized_counts.tsv",
        "results/deseq2/pca.pdf",
        "results/deseq2/ma_plot.pdf",
    ]
