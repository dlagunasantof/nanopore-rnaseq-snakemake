rule nanocount:
    """Transcript-level quantification of long reads with NanoCount (EM over multi-mapping reads)."""
    input:
        "results/alignments/{sample}.bam",
    output:
        "results/nanocount/{sample}.tsv",
    params:
        extra=config["nanocount"]["extra"],
    log:
        "logs/nanocount/{sample}.log",
    conda:
        "../envs/nanocount.yaml"
    shell:
        "NanoCount -i {input} -o {output} {params.extra} > {log} 2>&1"


rule merge_counts:
    """Merge per-sample tables into transcript and gene count matrices."""
    input:
        counts=expand("results/nanocount/{sample}.tsv", sample=SAMPLES),
        tx2gene="resources/reference/tx2gene.tsv",
    output:
        tx_counts="results/counts/transcript_counts.tsv",
        tx_tpm="results/counts/transcript_tpm.tsv",
        gene_counts="results/counts/gene_counts.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/merge_counts.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_counts.py"
