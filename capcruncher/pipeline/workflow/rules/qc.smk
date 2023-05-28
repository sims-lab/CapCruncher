import os
import pathlib


def get_fastq_basename(wildcards, output):
    return pathlib.Path(
        FASTQ_SAMPLES.translation[f"{wildcards.sample}_{wildcards.read}.fastq.gz"]
    ).stem.replace(".fastq", "")


rule fastqc:
    input:
        fq=lambda wc: FASTQ_SAMPLES.translation[f"{wc.sample}_{wc.read}.fastq.gz"],
    output:
        qc="capcruncher_output/qc/fastqc/{sample}_{read}_fastqc.html",
    params:
        outdir="capcruncher_output/qc/fastqc",
        tmpdir="capcruncher_output/qc/fastqc/{sample}_{read}",
        basename=lambda wc, output: get_fastq_basename(wc, output),
    threads: 12
    resources:
        mem_mb=1024,
    log:
        "capcruncher_output/logs/fastqc/{sample}_{read}.log",
    shell:
        """
        mkdir -p {params.tmpdir} &&
        fastqc -o {params.tmpdir} {input.fq} -t {threads} > {log} 2>&1 &&
        mv {params.tmpdir}/{params.basename}_fastqc.html {output.qc} &&
        mv {params.tmpdir}/{params.basename}_fastqc.zip {params.outdir}/{wildcards.sample}_{wildcards.read}_fastqc.zip &&
        rm -r {params.tmpdir}
        """


rule samtools_stats:
    input:
        bam="capcruncher_output/aligned/{sample}.bam",
        bai="capcruncher_output/aligned/{sample}.bam.bai",
    output:
        stats="capcruncher_output/qc/alignment_raw/{sample}.txt",
    threads: 1
    resources:
        mem_mb=1000,
    shell:
        """samtools stats {input.bam} > {output.stats}"""


rule multiqc:
    input:
        expand(
            "capcruncher_output/qc/fastqc/{sample}_{read}_fastqc.html",
            sample=SAMPLE_NAMES,
            read=[1, 2],
        ),
        expand("capcruncher_output/qc/alignment_raw/{sample}.txt", sample=SAMPLE_NAMES),
    output:
        "capcruncher_output/qc/full_qc_report.html",
    log:
        "capcruncher_output/logs/multiqc.log",
    resources:
        mem_mb=1000,
    shell:
        "multiqc -o capcruncher_output/qc capcruncher_output/qc -n full_qc_report.html --force > {log} 2>&1"
