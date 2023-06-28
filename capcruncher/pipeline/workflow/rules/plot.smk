
rule heatmaps_make_templates:
    input:
        counts=expand("capcruncher_output/{sample}/{sample}.hdf5", sample=SAMPLE_NAMES),
    output:
        "capcruncher_output/heatmap_templates/template_{viewpoint}_{bin_size}.yml",
    params:
        viewpoint="{viewpoint}",
        genes=config["plot"]["genes"] if os.path.exists(config["plot"]["genes"]) else "",
        bin_size="{bin_size}",
    shell:
        """
        capcruncher \
        plot \
        make-template \
        {input.counts} \
        {params.genes} \
        -o {output} \
        -v {params.viewpoint} \
        -b {params.bin_size}
        """


# rule pileups_make_template:
#     input:
#     output:
#         "capcruncher_output/pileup_templates/template_{viewpoint}.yml",
