rule minimap2_align:
    """
    Align reads to the transcriptome. Secondary alignments (-N) are kept on purpose:
    NanoCount uses them to assign multi-mapping reads to isoforms with an EM algorithm.
    The BAM is left unsorted, as NanoCount expects reads grouped by name.
    """
    input:
        reads="results/reads/{sample}.fastq.gz",
        index="resources/reference/transcriptome.mmi",
    output:
        "results/alignments/{sample}.bam",
    params:
        extra=config["minimap2"]["extra"],
    threads: config["threads"]["minimap2"]
    log:
        "logs/minimap2/{sample}.log",
    conda:
        "../envs/align.yaml"
    shell:
        "(minimap2 -t {threads} -ax map-ont {params.extra} {input.index} {input.reads} "
        "| samtools view -b -o {output} -) 2> {log}"


rule samtools_flagstat:
    """Mapping statistics, collected by MultiQC."""
    input:
        "results/alignments/{sample}.bam",
    output:
        "results/alignments/{sample}.flagstat.txt",
    log:
        "logs/flagstat/{sample}.log",
    conda:
        "../envs/align.yaml"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"
