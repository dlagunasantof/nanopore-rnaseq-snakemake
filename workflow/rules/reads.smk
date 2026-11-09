rule prepare_reads:
    """Copy reads into the results folder, optionally subsampling them (subsample_reads > 0)."""
    input:
        get_fastq,
    output:
        "results/reads/{sample}.fastq.gz",
    params:
        n=config.get("subsample_reads", 0),
        seed=config.get("seed", 42),
    log:
        "logs/prepare_reads/{sample}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        """
        if [ {params.n} -gt 0 ]; then
            seqkit sample -2 -s {params.seed} -n {params.n} {input} -o {output} 2> {log}
        else
            seqkit seq {input} -o {output} 2> {log}
        fi
        """


rule nanoplot:
    """Read-level QC (length, quality, yield) with NanoPlot."""
    input:
        "results/reads/{sample}.fastq.gz",
    output:
        stats="results/qc/nanoplot/{sample}/{sample}_NanoStats.txt",
        report="results/qc/nanoplot/{sample}/{sample}_NanoPlot-report.html",
    params:
        outdir=lambda wc, output: os.path.dirname(output.stats),
    threads: 2
    log:
        "logs/nanoplot/{sample}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        "NanoPlot -t {threads} --fastq {input} -o {params.outdir} "
        "--prefix {wildcards.sample}_ --tsv_stats --no_static > {log} 2>&1"
