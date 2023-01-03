from collections import defaultdict


def get_digestion_statistics(wc):
    stat_types = {
        "read_level_stats": "digestion.read.summary.csv",
        "histogram_unfiltered": "digestion.unfiltered.histogram.csv",
        "histogram_filtered": "digestion.filtered.histogram.csv",
    }

    stat_prefixes = []
    for sample in SAMPLE_NAMES:
        for part in get_fastq_partition_numbers_for_sample(wc, sample_name=sample):
            for combined in ["flashed", "pe"]:
                stat_prefixes.append(
                    f"capcruncher_statistics/digestion/data/{sample}_part{part}_{combined}."
                )

    stat_files = defaultdict(list)
    for stat_type, stat_suffix in stat_types.items():
        for stat_prefix in stat_prefixes:
            stat_files[stat_type].append(stat_prefix + stat_suffix)

    return stat_files


def get_filtering_statistics(wc):
    stat_types = {
        "read_level_stats": "read.stats.csv",
        "slice_level_stats": "slice.stats.csv",
    }

    stat_prefixes = []
    for sample in SAMPLE_NAMES:
        for part in get_fastq_partition_numbers_for_sample(wc, sample_name=sample):
            for combined in ["flashed", "pe"]:
                stat_prefixes.append(
                    f"capcruncher_statistics/filtering/data/{sample}_part{part}_{combined}."
                )

    stat_files = defaultdict(list)
    for stat_type, stat_suffix in stat_types.items():
        for stat_prefix in stat_prefixes:
            stat_files[stat_type].append(stat_prefix + stat_suffix)

    return stat_files


rule combine_stats_fastq_deduplication:
    input:
        fastq_deduplication=lambda wc: [
            f"capcruncher_statistics/deduplication/data/{sample}_part{part}.deduplication.csv"
            for sample in SAMPLE_NAMES
            for part in get_fastq_partition_numbers_for_sample(wc, sample_name=sample)
        ],
    output:
        "capcruncher_statistics/deduplication/fastq_deduplication.csv",
    script:
        "scripts/combine_deduplication_stats.py"


rule combine_stats_digestion:
    input:
        unpack(get_digestion_statistics),
    output:
        read_data="capcruncher_statistics/digestion/fastq_digestion.csv",
        histogram="capcruncher_statistics/digestion/fastq_digestion.histogram.csv",
    script:
        "scripts/combine_digestion_stats.py"


rule combine_stats_filtering:
    input:
        unpack(get_filtering_statistics),
    output:
        read_data="capcruncher_statistics/filtering/alignment_filtering.csv",
        slice_data="capcruncher_statistics/filtering/alignment_filtering_slice.csv",
    script:
        "scripts/combine_filtering_stats.py"


rule combine_stats_alignment_deduplication:
    input:
        read_level_stats=expand(
            "capcruncher_statistics/deduplication_by_coordinate/data/{sample}_{combined}.read.stats.csv",
            sample=SAMPLE_NAMES,
            combined=["flashed", "pe"],
        ),
    output:
        read_data="capcruncher_statistics/deduplication_by_coordinate/alignment_deduplication.csv",
    script:
        "scripts/combine_alignment_deduplication_stats.py"


rule merge_stats_filtering_and_alignment_deduplication:
    input:
        filtering=rules.combine_stats_filtering.output.read_data,
        alignment_deduplication=rules.combine_stats_alignment_deduplication.output.read_data,
    output:
        "capcruncher_statistics/filtering_and_alignment_deduplication.csv",
    log:
        "logs/merge_stats_filtering_and_alignment_deduplication.log",
    shell:
        """
        cat {input.filtering} > {output}
        cat {input.alignment_deduplication} | sed '1d' >> {output}
        """


rule combine_stats_cis_and_trans:
    input:
        cis_and_trans_stats=expand(
            "capcruncher_statistics/cis_and_trans_reporters/data/{sample}_{combined}.reporter.stats.csv",
            sample=SAMPLE_NAMES,
            combined=["flashed", "pe"],
        ),
    output:
        cis_and_trans_stats="capcruncher_statistics/cis_and_trans_reporters/cis_and_trans_reporters.csv",
    script:
        "scripts/combine_cis_and_trans_stats.py"


rule combine_stats_read_level:
    input:
        [
            rules.combine_stats_fastq_deduplication.output[0],
            rules.combine_stats_digestion.output.read_data,
            rules.merge_stats_filtering_and_alignment_deduplication.output[0],
        ],
    output:
        "capcruncher_statistics/run_statistics.csv",
    script:
        "scripts/combine_stats_read_level.py"


rule make_report:
    input:
        fastq_deduplication=rules.combine_stats_fastq_deduplication.output[0],
        digestion_read=rules.combine_stats_digestion.output.read_data,
        digestion_histogram=rules.combine_stats_digestion.output.histogram,
        reporters=rules.merge_stats_filtering_and_alignment_deduplication.output[0],
        cis_and_trans_stats=rules.combine_stats_cis_and_trans.output.cis_and_trans_stats,
        read_level_stats=rules.combine_stats_read_level.output[0],
    output:
        "capcruncher_statistics/capcruncher_report.html",
    script:
        "scripts/make_report.py"
