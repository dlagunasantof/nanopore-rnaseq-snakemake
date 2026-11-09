rule deseq2:
    """Gene-level differential expression between two conditions with DESeq2."""
    input:
        counts="results/counts/gene_counts.tsv",
        samples=config["samples"],
    output:
        results="results/deseq2/deseq2_results.tsv",
        normalized="results/deseq2/normalized_counts.tsv",
        pca="results/deseq2/pca.pdf",
        ma="results/deseq2/ma_plot.pdf",
    params:
        contrast=config["deseq2"]["contrast"],
        alpha=config["deseq2"]["alpha"],
        min_count=config["deseq2"]["min_count"],
    log:
        "logs/deseq2.log",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2.R"
