
rule regenerate_fastq:
    input:
        fq1=rules.fastq_rename.output.fq1,
        fq2=rules.fastq_rename.output.fq2,
        parquet="capcruncher_output/results/{sample}/{sample}.parquet",
    output:
        fq1="capcruncher_output/results/{sample}/{sample}_1.fastq.gz",
        fq2="capcruncher_output/results/{sample}/{sample}_2.fastq.gz",
    log:
        "capcruncher_output/logs/regenerate_fastq/{sample}.log",
    params:
        output_prefix=lambda wc, output: pathlib.Path(output.fq1).parent / wc.sample,
    shell:
        """
        capcruncher utilities regenerate-fastq \
            -p {input.parquet} \
            --output-prefix {params.output_prefix} \
            --fastq1 {input.fq1} \
            --fastq2 {input.fq2} \
            >{log} 2>&1
        """
