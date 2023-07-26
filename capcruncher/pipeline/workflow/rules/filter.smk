import capcruncher.pipeline.utils


def get_filtered_slices(wildcards):
    slices = dict()
    for combined_type in ["flashed", "pe"]:
        parts = get_rebalanced_parts(wildcards, combined=combined_type)
        slices[combined_type] = [
            f"capcruncher_output/interim/filtering/initial/{wildcards.sample}/{wildcards.sample}_part{part}_{combined_type}.slices.parquet"
            for part in parts
        ]
    return slices


rule filter_alignments:
    input:
        bam=rules.align_bowtie2.output.bam,
        annotations=rules.annotate.output.annotated,
    output:
        filtered_slices="capcruncher_output/interim/filtering/initial/{sample}/{sample}_part{part}_{combined}.slices.parquet",
        stats_read="capcruncher_output/interim/statistics/filtering/data/{sample}_part{part}_{combined}.read.stats.csv",
        stats_slice="capcruncher_output/interim/statistics/filtering/data/{sample}_part{part}_{combined}.slice.stats.csv",
    params:
        analysis_method=config["analysis"]["method"],
        sample_name=lambda wildcards, output: wildcards.sample,
        output_prefix=lambda wildcards, output: output.filtered_slices.replace(
            ".slices.parquet", ""
        ),
        stats_prefix=lambda wildcards, output: output.stats_read.replace(
            ".read.stats.csv", ""
        ),
        read_type=lambda wildcards, output: wildcards.combined,
        custom_filtering=capcruncher.pipeline.utils.validate_custom_filtering(config),
    resources:
        mem_mb=5000,
    log:
        "capcruncher_output/interim/statistics/filtering/logs/{sample}_part{part}_{combined}.log",
    shell:
        """
        capcruncher \
        alignments \
        filter \
        {params.analysis_method} \
        -b {input.bam} \
        -a {input.annotations} \
        -o {params.output_prefix} \
        --stats-prefix {params.stats_prefix} \
        --sample-name {params.sample_name} \
        --read-type {params.read_type} \
        --no-cis-and-trans-stats \
        --no-fragments \
        {params.custom_filtering} > {log} 2>&1
        """


rule split_flashed_and_pe_datasets:
    input:
        unpack(get_filtered_slices),
    output:
        slices_flashed=directory(
            "capcruncher_output/interim/filtering/repartitioned/{sample}/flashed/"
        ),
        slices_pe=directory(
            "capcruncher_output/interim/filtering/repartitioned/{sample}/pe/"
        ),
    shell:
        """
        mkdir -p {output.slices_flashed}
        mkdir -p {output.slices_pe}
        mv {input.flashed} {output.slices_flashed}
        mv {input.pe} {output.slices_pe}
        """


rule remove_duplicate_coordinates:
    input:
        slices_directory="capcruncher_output/interim/filtering/repartitioned/{sample}/{combined}/",
    output:
        slices=directory(
            "capcruncher_output/interim/filtering/deduplicated/{sample}/{combined}"
        ),
        stats_read="capcruncher_output/interim/statistics/deduplication_by_coordinate/data/{sample}_{combined}.read.stats.csv",
    params:
        sample_name=lambda wildcards, output: wildcards.sample,
        stats_prefix=lambda wildcards, output: output.stats_read.replace(
            ".read.stats.csv", ""
        ),
        read_type=lambda wildcards, output: wildcards.combined,
    resources:
        mem_mb=lambda wc, attempt: 3000 * 2**attempt,
    threads: 12
    log:
        "capcruncher_output/logs/remove_duplicate_coordinates/{sample}_{combined}.log",
    script:
        "../scripts/remove_duplicate_coordinates.py"


rule combine_flashed_and_pe_post_deduplication:
    input:
        slices=expand(
            "capcruncher_output/interim/filtering/deduplicated/{{sample}}/{combined}",
            combined=["flashed", "pe"],
        ),
    output:
        slices=directory("capcruncher_output/results/{sample}/{sample}.parquet"),
    shell:
        """
        mkdir -p {output.slices}
        mv {input.slices[0]}/*.parquet {output.slices} || mkdir -p {output.slices}
        mv {input.slices[1]}/*.parquet {output.slices} || mkdir -p {output.slices}
        """


rule cis_and_trans_stats:
    input:
        slices="capcruncher_output/results/{sample}/{sample}.parquet",
    output:
        stats="capcruncher_output/interim/statistics/cis_and_trans_reporters/data/{sample}.reporter.stats.csv",
    params:
        sample_name=lambda wildcards, output: wildcards.sample,
        analysis_method=config["analysis"]["method"],
    resources:
        mem_mb=3000,
    log:
        "capcruncher_output/logs/cis_and_trans_stats/{sample}.log",
    shell:
        """
        capcruncher \
        utilities \
        cis-and-trans-stats \
        {input.slices} \
        --sample-name {params.sample_name} \
        -o {output.stats} \
        """
