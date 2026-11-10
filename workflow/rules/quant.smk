rule oarfish:
    """
    Transcript-level quantification of long reads with oarfish (EM over multi-mapping reads).
    Alignments on both strands are kept: PCR-cDNA reads come in both orientations, unlike
    direct RNA reads, so a strand filter would discard real reads.
    """
    input:
        "results/alignments/{sample}.bam",
    output:
        quant="results/oarfish/{sample}.quant",
    params:
        prefix=lambda wc, output: output.quant.removesuffix(".quant"),
        extra=config["oarfish"]["extra"],
    threads: 4
    log:
        "logs/oarfish/{sample}.log",
    conda:
        "../envs/oarfish.yaml"
    shell:
        "oarfish -j {threads} -a {input} -o {params.prefix} {params.extra} > {log} 2>&1"


rule merge_counts:
    """Merge per-sample tables into transcript and gene count matrices."""
    input:
        counts=expand("results/oarfish/{sample}.quant", sample=SAMPLES),
        tx2gene="resources/reference/tx2gene.tsv",
    output:
        tx_counts="results/counts/transcript_counts.tsv",
        tx_cpm="results/counts/transcript_cpm.tsv",
        gene_counts="results/counts/gene_counts.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/merge_counts.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_counts.py"
