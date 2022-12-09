import capcruncher.pipeline.utils as pipeline_utils


rule fastqc:
    input:
        fq = "FASTQ",
    output:
        qc = "OUTPUT",
    params:
        outdir = "OUTDIR",
    threads:
        4,
    shell:
        """
        fastqc -q -t {threads} --nogroup --outdir {params.outdir} {input.fq}
        """

rule multiqc:
    input:
        reports = "INPUTS",
    output:
        report = "OUTPUT",
    threads:
        1,
    log:
        "LOG"
    run:

        outdir = os.path.dirname(output.report)
        basename = os.path.basename(output.report)

        if not isinstance(input.reports, list):
            reports = [*input.reports]
        else:
            reports = input
        
        dirnames = [os.path.dirname(x) for x in reports]
        dirnames = list(set(dirnames))
        search_dirs = " ".join(dirnames)

        cmd = f"multiqc -n {basename} -o {outdir} {search_dirs} --force > {log} 2>&1"

        if workflow.use_singularity:
            cmd = pipeline_utils.get_singularity_command(command=cmd,
                                                workflow=workflow,)
        shell(cmd)


use rule fastqc as fastqc_raw with:
    input:
        fq = "{sample}.fastq.gz",
    output:
        qc = "qc/fastq_raw/{sample}_fastqc.html",
    params:
        outdir = "qc/fastq_raw"

use rule fastqc as fastqc_trimmed with:
    input:
        fq = "{sample}.fastq.gz",
    output:
        qc = "qc/fastq_trimmed/{sample}_fastqc.html",
    params:
        outdir = "qc/fastq_trimmed"

use rule multiqc as multiqc_fastq_raw with:
    input:
        reports = expand("qc/fastq_raw/{sample}_{read}_fastqc.html", sample=SAMPLE_NAMES, read=[1,2]),
    output:
        report = "qc/fastq_qc_raw_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_raw_report.log"

use rule multiqc as multiqc_fastq_trimmed with:
    input:
        reports = expand("qc/fastq_trimmed/{sample}_{read}_fastqc.html", sample=SAMPLE_NAMES, read=[1,2]),
    output:
        report = "qc/fastq_qc_trimmed_report.html"
    threads:
        1,
    log:
        "logs/qc/fastq_qc_trimmed_report.log"