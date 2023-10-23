import os
import pathlib
import capcruncher.pipeline.utils


rule fastqc:
    input:
        fq="capcruncher_output/interim/fastq/{sample}_{read}.fastq.gz",
    output:
        qc="capcruncher_output/interim/qc/fastqc/{sample}_{read}_fastqc.html",
    params:
        outdir=lambda wc, output: str(pathlib.Path(output.qc).parent),
        tmpdir="capcruncher_output/interim/qc/fastqc/{sample}_{read}",
    threads: 1
    resources:
        mem_mb=1024,
    log:
        "capcruncher_output/logs/fastqc/{sample}_{read}.log",
    shell:
        """
        mkdir -p {params.tmpdir} &&
        fastqc -o {params.tmpdir} {input.fq} -t {threads} > {log} 2>&1 &&
        mv {params.tmpdir}/{wildcards.sample}_{wildcards.read}_fastqc.html {output.qc} &&
        mv {params.tmpdir}/{wildcards.sample}_{wildcards.read}_fastqc.zip {params.outdir}/{wildcards.sample}_{wildcards.read}_fastqc.zip &&
        rm -r {params.tmpdir}
        """


rule samtools_stats:
    input:
        bam="capcruncher_output/results/{sample}/{sample}.bam",
        bai="capcruncher_output/results/{sample}/{sample}.bam.bai",
    output:
        stats=temp("capcruncher_output/interim/qc/alignment_raw/{sample}.txt"),
    threads: 1
    resources:
        mem_mb=1000,
    shell:
        """samtools stats {input.bam} > {output.stats}"""


rule multiqc_report:
    input:
        expand(
            "capcruncher_output/interim/qc/fastqc/{sample}_{read}_fastqc.html",
            sample=SAMPLE_NAMES,
            read=[1, 2],
        ),
        expand(
            "capcruncher_output/interim/qc/alignment_raw/{sample}.txt",
            sample=SAMPLE_NAMES,
        ),
    output:
        "capcruncher_output/results/full_qc_report.html",
    log:
        "capcruncher_output/logs/multiqc.log",
    params:
        outdir=lambda wc, output: str(pathlib.Path(output[0]).parent),
        dir_analysis="capcruncher_output/interim/qc",
    resources:
        mem_mb=1000,
    shell:
        "multiqc -o {params.outdir} {params.dir_analysis} -n full_qc_report.html --force > {log} 2>&1"


rule multiqc_full:
    input:
        report="capcruncher_output/results/full_qc_report.html",
    output:
        report="capcruncher_output/interim/statistics/multiqc_full_data/multiqc_report.html",
        trimming_data="capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_cutadapt.txt",
        flash_data="capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_flash_combo_stats.txt",
        bowtie2_data="capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_bowtie2.txt",
    log:
        "capcruncher_output/logs/multiqc_full.log",
    shell:
        "multiqc capcruncher_output/ --outdir capcruncher_output/interim/statistics/multiqc_full_data --force > {log} 2>&1"
