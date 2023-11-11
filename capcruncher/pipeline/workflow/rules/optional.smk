
rule regenerate_fastq:
    input:
        fq1=rules.fastq_rename.output.fq1,
        fq2=rules.fastq_rename.output.fq2,
        parquet="capcruncher_output/results/{sample}/{sample}.parquet"
    output:
        fq1="capcruncher_output/results/{sample}/{sample}_1.fastq.gz",
        fq2="capcruncher_output/results/{sample}/{sample}_2.fastq.gz",
    params:
        output_prefix="capcruncher_output/results/{sample}/{sample}",
    shell:
        """
        capcruncher regenerate-fastq \
            -p {input.parquet} \
            --output-prefix {params.output_prefix} \
            --fastq1 {input.fq1} \
            --fastq2 {input.fq2}
        """
        