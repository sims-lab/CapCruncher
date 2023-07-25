from collections import defaultdict
import capcruncher.pipeline.utils


def get_digestion_statistics(wc, sample_names: List[str]):
    stat_types = {
        "read_level_stats": "digestion.read.summary.csv",
        "histogram_unfiltered": "digestion.unfiltered.histogram.csv",
        "histogram_filtered": "digestion.filtered.histogram.csv",
    }

    stat_prefixes = []
    for sample in sample_names:
        for combined in ["flashed", "pe"]:
            for part in get_rebalanced_parts(wc, combined=combined, sample=sample):
                stat_prefixes.append(
                    f"capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_{combined}."
                )

    stat_files = defaultdict(list)
    for stat_type, stat_suffix in stat_types.items():
        for stat_prefix in stat_prefixes:
            stat_files[stat_type].append(stat_prefix + stat_suffix)

    return stat_files


def get_filtering_statistics(wc, sample_names: List[str]):
    stat_types = {
        "read_level_stats": "read.stats.csv",
        "slice_level_stats": "slice.stats.csv",
    }

    stat_prefixes = []
    for sample in sample_names:
        for combined in ["flashed", "pe"]:
            for part in get_rebalanced_parts(wc, combined=combined, sample=sample):
                stat_prefixes.append(
                    f"capcruncher_output/interim/statistics/filtering/data/{sample}_part{part}_{combined}."
                )

    stat_files = defaultdict(list)
    for stat_type, stat_suffix in stat_types.items():
        for stat_prefix in stat_prefixes:
            stat_files[stat_type].append(stat_prefix + stat_suffix)

    return stat_files


def get_stat_parts(wc, sample_names: List[str]):
    files = []
    for sample in sample_names:
        for part in get_fastq_split_1(wc):
            files.append(
                f"capcruncher_output/interim/statistics/deduplication/data/{sample}_part{part}.deduplication.csv"
            )
    return files


if not CAPCRUNCHER_TOOLS:

    rule combine_stats_fastq_deduplication:
        input:
            fastq_deduplication=get_stat_parts,
        output:
            "capcruncher_output/interim/statistics/deduplication/fastq_deduplication.csv",
        script:
            "../scripts/combine_deduplication_stats.py"

else:

    rule combine_stats_fastq_deduplication:
        input:
            fastq_deduplication=expand(
                "capcruncher_output/interim/statistics/deduplication/data/{sample}.deduplication.csv",
                sample=SAMPLE_NAMES,
            ),
        output:
            "capcruncher_output/interim/statistics/deduplication/fastq_deduplication.csv",
        script:
            "../scripts/combine_deduplication_stats.py"


rule combine_stats_digestion:
    input:
        unpack(lambda wc: get_digestion_statistics(wc)),
    output:
        read_data="capcruncher_output/interim/statistics/digestion/fastq_digestion.csv",
        histogram="capcruncher_output/interim/statistics/digestion/fastq_digestion.histogram.csv",
    script:
        "../scripts/combine_digestion_stats.py"


rule combine_stats_filtering:
    input:
        unpack(get_filtering_statistics),
    output:
        read_data="capcruncher_output/interim/statistics/filtering/alignment_filtering.csv",
        slice_data="capcruncher_output/interim/statistics/filtering/alignment_filtering_slice.csv",
    script:
        "../scripts/combine_filtering_stats.py"


rule combine_stats_alignment_deduplication:
    input:
        read_level_stats=expand(
            "capcruncher_output/interim/statistics/deduplication_by_coordinate/data/{sample}_{combined}.read.stats.csv",
            sample=SAMPLE_NAMES,
            combined=["flashed", "pe"],
        ),
    output:
        read_data="capcruncher_output/interim/statistics/deduplication_by_coordinate/alignment_deduplication.csv",
    script:
        "../scripts/combine_alignment_deduplication_stats.py"


rule merge_stats_filtering_and_alignment_deduplication:
    input:
        filtering=rules.combine_stats_filtering.output.read_data,
        alignment_deduplication=rules.combine_stats_alignment_deduplication.output.read_data,
    output:
        "capcruncher_output/interim/statistics/filtering_and_alignment_deduplication.csv",
    log:
        "capcruncher_output/logs/merge_stats_filtering_and_alignment_deduplication.log",
    shell:
        """
        cat {input.filtering} > {output}
        cat {input.alignment_deduplication} | sed '1d' >> {output}
        """


rule combine_stats_cis_and_trans:
    input:
        cis_and_trans_stats=expand(
            "capcruncher_output/interim/statistics/cis_and_trans_reporters/data/{sample}.reporter.stats.csv",
            sample=SAMPLE_NAMES,
        ),
    output:
        cis_and_trans_stats="capcruncher_output/interim/statistics/cis_and_trans_reporters/cis_and_trans_reporters.csv",
    script:
        "../scripts/combine_cis_and_trans_stats.py"


rule combine_stats_read_level:
    input:
        [
            rules.combine_stats_fastq_deduplication.output[0],
            rules.combine_stats_digestion.output.read_data,
            rules.merge_stats_filtering_and_alignment_deduplication.output[0],
        ],
    output:
        "capcruncher_output/interim/statistics/run_statistics.csv",
    script:
        "../scripts/combine_stats_read_level.py"


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


rule make_report:
    input:
        template=rules.copy_report_template.output[0],
        fastq_deduplication=rules.combine_stats_fastq_deduplication.output[0],
        digestion_read=rules.combine_stats_digestion.output.read_data,
        digestion_histogram=rules.combine_stats_digestion.output.histogram,
        reporters=rules.merge_stats_filtering_and_alignment_deduplication.output[0],
        cis_and_trans_stats=rules.combine_stats_cis_and_trans.output.cis_and_trans_stats,
        read_level_stats=rules.combine_stats_read_level.output[0],
    output:
        "capcruncher_output/results/capcruncher_report.html",
    params:
        outdir=lambda wildcards, output: pathlib.Path(output[0]).parent,
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
        -P fastq_deduplication_path:$(realpath {input.fastq_deduplication}) \
        -P fastq_digestion_read_path:$(realpath {input.digestion_read}) \
        -P fastq_digestion_hist_path:$(realpath {input.digestion_histogram}) \
        -P reporter_read_path:$(realpath {input.reporters}) \
        -P reporter_cis_trans_path:$(realpath {input.cis_and_trans_stats}) \
        -P run_stats_path:$(realpath {input.read_level_stats}) \
        --log {log}

        rm {params.outdir}/capcruncher_report.qmd
        """
