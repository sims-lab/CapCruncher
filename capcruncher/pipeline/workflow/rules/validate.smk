def get_annotated_slices(wildcards):
    slices = dict()
    for combined_type in ["flashed", "pe"]:
        parts = get_rebalanced_parts(wildcards, combined=combined_type)
        slices[combined_type] = [
            f"capcruncher_output/interim/annotate/{wildcards.sample}/{wildcards.sample}_part{part}_{combined_type}.slices.parquet"
            for part in parts
        ]
    return [*slices["flashed"], *slices["pe"]]


rule check_n_bins_per_viewpoint:
    input:
        bins=rules.digest_genome.output.bed,
        viewpoints=config["analysis"]["viewpoints"],
    output:
        sentinel="capcruncher_output/resources/validation/check_n_bins_per_viewpoint.sentinel",
        n_bins_per_viewpoint="capcruncher_output/resources/validation/n_bins_per_viewpoint.tsv",
    params:
        ignore_multiple_bins_per_viewpoint=IGNORE_MULTIPLE_FRAGMENTS_PER_VIEWPOINT,
    script:
        "scripts/validation_check_n_bins_per_viewpoint.py"


rule confirm_annotated_viewpoints_present:
    input:
        slices=get_annotated_slices,
        viewpoints=config["analysis"]["viewpoints"],
    output:
        sentinel="capcruncher_output/resources/validation/confirm_annotated_viewpoints_present.sentinel",
        viewpoints_present="capcruncher_output/resources/validation/annotated_viewpoints_present.tsv",
    script:
        "scripts/validation_confirm_annotated_viewpoints_present.py"
