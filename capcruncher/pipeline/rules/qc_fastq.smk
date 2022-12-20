rule fastqc_raw:
    input:
        fq="{sample}_{read}.fastq.gz",
    output:
        qc="capcruncher_statistics/fastqc_raw/{sample}/{sample}_{read}_fastqc.html",
    params:
        outdir=lambda wc, output: os.path.dirname(output.qc),
    threads: 4
    shell:
        """
        fastqc -q -t {threads} --nogroup --outdir {params.outdir} {input.fq}
        """


rule multiqc_raw:
    input:
        fastqc=expand(
            "capcruncher_statistics/fastqc_raw/{sample}/{sample}_{read}_fastqc.html",
            sample=SAMPLE_NAMES,
            read=[1, 2],
        ),
    output:
        report="capcruncher_statistics/fastqc_raw/multiqc/fastqc_raw.html",
    params:
        input_dir="capcruncher_statistics/fastqc_raw/",
        output_dir="capcruncher_statistics/fastqc_raw/multiqc/",
        report_name="fastqc_raw.html",
    threads: 1
    log:
        "logs/qc/fastq_qc_raw_report.log",
    shell:
        """
        multiqc -n {params.report_name} -o {params.output_dir} {params.input_dir} --force > {log} 2>&1
        """
