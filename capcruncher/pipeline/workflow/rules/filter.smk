import capcruncher.pipeline.utils

# rule check_viewpoints_annotated:
#     input:
#         slices=get_annotated_slices,
#         viewpoints=config["analysis"]["viewpoints"],
#     output:
#         sentinel="capcruncher_output/resources/validation/{sample}.check_viewpoints.sentinel",
#         viewpoints_present="capcruncher_output/resources/validation/{sample}.annotated_viewpoints_present.tsv",
#     script:
#         "../scripts/validation_confirm_annotated_viewpoints_present.py"


rule filter_alignments:
    input:
        bam=rules.align_bowtie2.output.bam,
        annotations=rules.annotate.output.annotated,
    output:
        filtered_slices=temp(
            "capcruncher_output/interim/filtering/initial/{sample}/{sample}_part{part}_{combined}.slices.parquet"
        ),
        statistics="capcruncher_output/interim/statistics/filtering/data/{sample}_part{part}_{combined}.json",
    log:
        "capcruncher_output/interim/statistics/filtering/logs/{sample}_part{part}_{combined}.log",
    resources:
        mem=lambda wildcards, attempt: scale_memory(5, attempt),
    params:
        analysis_method=config["analysis"]["method"],
        sample_name=lambda wildcards, output: wildcards.sample,
        output_prefix=lambda wildcards, output: output.filtered_slices.replace(
            ".slices.parquet", ""
        ),
        read_type=lambda wildcards, output: wildcards.combined,
        filter_profile=capcruncher.pipeline.utils.validate_filter_profile(config),
    shell:
        """
        capcruncher \
            alignments \
            filter \
            {params.analysis_method} \
            -b {input.bam} \
            -a {input.annotations} \
            -o {params.output_prefix} \
            --statistics {output.statistics} \
            --sample-name {params.sample_name} \
            --read-type {params.read_type} \
            --no-fragments \
            {params.filter_profile} >{log} 2>&1
        """


rule split_flashed_and_pe_datasets:
    input:
        flashed=lambda wc: get_filtered_slices(wc)["flashed"],
        pe=lambda wc: get_filtered_slices(wc)["pe"],
    output:
        slices_flashed=temp(
            directory(
                "capcruncher_output/interim/filtering/repartitioned/{sample}/flashed/"
            )
        ),
        slices_pe=temp(
            directory(
                "capcruncher_output/interim/filtering/repartitioned/{sample}/pe/"
            )
        ),
    log:
        "capcruncher_output/logs/split_flashed_and_pe_datasets/{sample}.log",
    script:
        "../scripts/repartition_filtered_slices.py"


rule remove_duplicate_coordinates:
    input:
        slices_directory="capcruncher_output/interim/filtering/repartitioned/{sample}/{combined}/",
    output:
        slices=temp(
            directory(
                "capcruncher_output/interim/filtering/deduplicated/{sample}/{combined}"
            )
        ),
        statistics="capcruncher_output/interim/statistics/deduplication_final/data/{sample}_{combined}.json",
    log:
        "capcruncher_output/logs/remove_duplicate_coordinates/{sample}_{combined}.log",
    threads: 12
    resources:
        mem=lambda wildcards, attempt: scale_memory(3, attempt),
    params:
        sample_name=lambda wildcards, output: wildcards.sample,
        read_type=lambda wildcards, output: wildcards.combined,
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
    log:
        "capcruncher_output/logs/combine_flashed_and_pe_post_deduplication/{sample}.log",
    params:
        source_dir=lambda wc, input: pathlib.Path(input.slices[0]).parent,
        dest_dir=lambda wc, output: pathlib.Path(output.slices),
    script:
        "../scripts/combine_deduplicated_slices.py"


rule cis_and_trans_stats:
    input:
        slices="capcruncher_output/results/{sample}/{sample}.parquet",
    output:
        stats="capcruncher_output/interim/statistics/cis_and_trans_reporters/data/{sample}.json",
    log:
        "capcruncher_output/logs/cis_and_trans_stats/{sample}.log",
    resources:
        mem=lambda wildcards, attempt: scale_memory(3, attempt),
    params:
        sample_name=lambda wildcards, output: wildcards.sample,
        analysis_method=config["analysis"]["method"],
    shell:
        """
        capcruncher \
            utilities \
            cis-and-trans-stats \
            {input.slices} \
            --assay {params.analysis_method} \
            --sample-name {params.sample_name} \
            -o {output.stats} \
            >{log} 2>&1
        """


localrules:
    split_flashed_and_pe_datasets,
    combine_flashed_and_pe_post_deduplication,
