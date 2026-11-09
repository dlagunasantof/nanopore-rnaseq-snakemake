rule multiqc:
    """One HTML report with read QC (NanoPlot) and mapping statistics (samtools)."""
    input:
        expand("results/qc/nanoplot/{sample}/{sample}_NanoStats.txt", sample=SAMPLES),
        expand("results/alignments/{sample}.flagstat.txt", sample=SAMPLES),
    output:
        "results/multiqc/multiqc_report.html",
    params:
        outdir=lambda wc, output: os.path.dirname(output[0]),
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force --outdir {params.outdir} {input} > {log} 2>&1"
