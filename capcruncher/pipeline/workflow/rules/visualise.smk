import capcruncher.pipeline.utils


rule viewpoints_to_bigbed:
    input:
        viewpoints=config["analysis"]["viewpoints"],
    params:
        chrom_sizes=config["genome"]["chrom_sizes"],
    output:
        "capcruncher_output/resources/viewpoints/viewpoints.bigBed",
    log:
        "capcruncher_output/logs/viewpoints_to_bigbed.log",
    shell:
        """
        sort -k1,1 -k2,2n {input.viewpoints} > {output}.tmp
        bedToBigBed {output}.tmp {params.chrom_sizes} {output} > {log} 2>&1
        rm {output}.tmp
        """


rule create_ucsc_hub:
    input:
        viewpoints=rules.viewpoints_to_bigbed.output[0],
        bigwigs=expand(
            "capcruncher_output/results/{sample}/bigwigs/{norm}/{sample}_{viewpoint}.bigWig",
            sample=SAMPLE_NAMES,
            norm=["raw", "norm"],
            viewpoint=VIEWPOINT_NAMES,
        ),
        bigwigs_summary=expand(
            "capcruncher_output/results/comparisons/bigwigs/{group}.{method}-summary.{viewpoint}.bigWig",
            group=DESIGN["condition"].unique(),
            method=SUMMARY_METHODS,
            viewpoint=VIEWPOINT_NAMES,
        )
        if AGGREGATE_SAMPLES
        else [],
        bigwigs_comparison=expand(
            "capcruncher_output/results/comparisons/bigwigs/{comparison}.{method}-subtraction.{viewpoint}.bigWig",
            comparison=[
            f"{a}-{b}"
                for a, b in itertools.permutations(DESIGN["condition"].unique(), 2)
            ],
            method=SUMMARY_METHODS,
            viewpoint=VIEWPOINT_NAMES,
        )
        if COMPARE_SAMPLES
        else [],
        report=rules.make_report.output[0],
    output:
        directory(config["hub"]["dir"]),
    wildcard_constraints:
        comparison=r"[A-Za-z0-9_.]+-[A-Za-z0-9_.]+",
        group=r"[A-Za-z0-9_.]+",
    params:
        color_by=config["hub"].get("color_by", "sample"),
        genome=config["genome"]["name"],
        custom_genome=config["hub"].get("custom_genome", None),
        genome_twobit=config["genome"].get("twobit", None),
        hub_name=config["hub"].get("name"),
        hub_email=config["hub"].get("email"),
        genome_organism=config["genome"].get("organism"),
        genome_default_position=config["genome"].get("genome_default_position"),
    log:
        "capcruncher_output/logs/create_ucsc_hub.log",
    script:
        "../scripts/make_ucsc_hub.py"


rule plot:
    input:
        bigwigs=lambda wc: capcruncher.pipeline.utils.get_files_to_plot(
            wc, DESIGN, ASSAY, SAMPLE_NAMES, SUMMARY_METHODS, COMPARE_SAMPLES
        )["bigwigs"],
        subtractions=lambda wc: capcruncher.pipeline.utils.get_files_to_plot(
            wc, DESIGN, ASSAY, SAMPLE_NAMES, SUMMARY_METHODS, COMPARE_SAMPLES
        )["subtractions"],
        bigwigs_collection=lambda wc: capcruncher.pipeline.utils.get_files_to_plot(
            wc, DESIGN, ASSAY, SAMPLE_NAMES, SUMMARY_METHODS, COMPARE_SAMPLES
        )["bigwigs_collection"],
        heatmaps=lambda wc: capcruncher.pipeline.utils.get_files_to_plot(
            wc, DESIGN, ASSAY, SAMPLE_NAMES, SUMMARY_METHODS, COMPARE_SAMPLES
        )["heatmaps"],
        viewpoints=config["analysis"]["viewpoints"],
    output:
        template="capcruncher_output/results/figures/{viewpoint}.toml",
        fig="capcruncher_output/results/figures/{viewpoint}.pdf",
    params:
        coordinates=lambda wc: capcruncher.pipeline.utils.get_plotting_coordinates(
            wc, config
        ),
        viewpoint="{viewpoint}",
        design=DESIGN,
        genes=config["plot"].get("genes", ""),
        binsize=config["analysis"].get("bin_sizes", [None])[0],
        normalization_method=config["plot"].get("normalisation", "raw"),
    wildcard_constraints:
        comparison=r"[A-Za-z0-9_.]+-[A-Za-z0-9_.]+",
        group=r"[A-Za-z0-9_.]+",
    log:
        "capcruncher_output/logs/plot/{viewpoint}.log",
    threads: 1
    script:
        "../scripts/plot.py"


localrules:
    create_ucsc_hub,
    plot,
