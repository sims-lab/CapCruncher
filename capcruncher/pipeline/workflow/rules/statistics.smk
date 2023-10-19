from collections import defaultdict
import capcruncher.pipeline.utils
from typing import List


def get_digestion_statistics(wc, sample_names: List[str]):
    
    stat_files = []
    for sample in sample_names:
        for combined in ["flashed", "pe"]:
            for part in get_rebalanced_parts(wc, combined=combined, sample=sample):
                stat_files.append(
                    f"capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_{combined}.json"
                )
    
    return stat_files

def get_filtering_statistics(wc, sample_names: List[str]):

    stat_files = []
    for sample in sample_names:
        for combined in ["flashed", "pe"]:
            for part in get_rebalanced_parts(wc, combined=combined, sample=sample):
                stat_files.append(
                    f"capcruncher_output/interim/statistics/filtering/data/{sample}_part{part}_{combined}.json"
                )
    
    return stat_files


rule copy_report_template:
    input:
        template=workflow.source_path("../report/capcruncher_report.qmd"),
    output:
        "capcruncher_output/results/capcruncher_report.qmd",
    container:
        None
    shell:
        """
        cp {input.template} {output}
        """

rule extract_trimming_data:
    input:
        rules.multiqc_full.output.trimming_data,
    output:
        "capcruncher_output/interim/statistics/trimming/trimming.json",
    script:
        "../scripts/extract_trimming_data.py"

rule extract_flash_data:
    input:
        rules.multiqc_full.output.flash_data,
    output:
        "capcruncher_output/interim/statistics/flash/flash.json",
    script:
        "../scripts/extract_flash_data.py"


rule make_report:
    input:
        template=rules.copy_report_template.output[0],
        fastq_deduplication=expand(
            "capcruncher_output/interim/statistics/deduplication/data/{sample}.deduplication.json",
            sample=SAMPLE_NAMES,
        ),
        fastq_trimming=rules.extract_trimming_data.output[0],
        fastq_flash=rules.extract_flash_data.output[0],
        fastq_digestion=lambda wc: get_digestion_statistics(wc, SAMPLE_NAMES),
        reporters=lambda wc: get_filtering_statistics(wc, SAMPLE_NAMES),
        cis_and_trans_stats=expand(
            "capcruncher_output/interim/statistics/cis_and_trans_reporters/data/{sample}.json",
            sample=SAMPLE_NAMES,
        ),
    output:
        "capcruncher_output/results/capcruncher_report.html",
    params:
        outdir=lambda wildcards, output: pathlib.Path(output[0]).parent,
        fastq_deduplication_path="capcruncher_output/interim/statistics/deduplication/data/",
        fastq_digestion_path="capcruncher_output/interim/statistics/digestion/data/",
        reporter_filtering_path="capcruncher_output/interim/statistics/filtering/data/",
        reporter_cis_trans_path="capcruncher_output/interim/statistics/cis_and_trans_reporters/data/",
    log:
        "capcruncher_output/logs/make_report.log",
    shell:
        """
        export XDG_RUNTIME_DIR=$(mktemp -d);
        quarto \
        render \
        {params.outdir}/capcruncher_report.qmd \
        --to html \
        --execute \
        -P fastq_deduplication_path:$(realpath {params.fastq_deduplication_path}) \
        -P fastq_trimming_path:$(realpath {input.fastq_trimming}) \
        -P fastq_flash_path:$(realpath {input.fastq_flash}) \
        -P fastq_digestion_path:$(realpath {params.fastq_digestion_path}) \
        -P reporter_filtering_path:$(realpath {params.reporter_filtering_path}) \
        -P reporter_cis_trans_path:$(realpath {params.reporter_cis_trans_path}) \
        --log {log} \
        2> {log}.err;

        rm {params.outdir}/capcruncher_report.qmd
        """
    
