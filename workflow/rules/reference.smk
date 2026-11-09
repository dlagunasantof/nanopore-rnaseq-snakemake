rule prepare_reference:
    """Clean transcript IDs in the FASTA and build a transcript-to-gene table from the GTF."""
    input:
        fasta=config["reference"]["transcriptome_fasta"],
        gtf=config["reference"]["annotation_gtf"],
    output:
        fasta="resources/reference/transcriptome.fa",
        tx2gene="resources/reference/tx2gene.tsv",
    log:
        "logs/reference/prepare_reference.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/prepare_reference.py"


rule minimap2_index:
    """Index the transcriptome once so every sample reuses it."""
    input:
        "resources/reference/transcriptome.fa",
    output:
        "resources/reference/transcriptome.mmi",
    threads: config["threads"]["minimap2"]
    log:
        "logs/reference/minimap2_index.log",
    conda:
        "../envs/align.yaml"
    shell:
        "minimap2 -t {threads} -x map-ont -d {output} {input} > {log} 2>&1"
