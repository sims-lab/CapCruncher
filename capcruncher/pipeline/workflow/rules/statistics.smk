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
        fastq_deduplication=expand(
            "capcruncher_output/interim/statistics/deduplication/data/{sample}.deduplication.json",
            sample=SAMPLE_NAMES,
        ),
        fastq_trimming=rules.extract_trimming_data.output[0],
        fastq_flash=rules.extract_flash_data.output[0],
        fastq_digestion=lambda wc: get_digestion_statistics(wc, SAMPLE_NAMES),
        reporters_filtering=lambda wc: get_filtering_statistics(wc, SAMPLE_NAMES),
        reporters_deduplication=expand(
            "capcruncher_output/interim/statistics/deduplication_final/data/{sample}_{combined}.json",
            sample=SAMPLE_NAMES,
            combined=["flashed", "pe"],
        ),
        cis_and_trans_stats=expand(
            "capcruncher_output/interim/statistics/cis_and_trans_reporters/data/{sample}.json",
            sample=SAMPLE_NAMES,
        ),
    output:
        "capcruncher_output/results/capcruncher_report.html",
    log:
        "capcruncher_output/logs/make_report.log",
    script:
        "../report/make_report.py"
