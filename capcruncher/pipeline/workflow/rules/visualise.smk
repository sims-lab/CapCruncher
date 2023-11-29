import capcruncher.pipeline.utils


rule viewpoints_to_bigbed:
    input:
        viewpoints=config["analysis"]["viewpoints"],
    params:
        chrom_sizes=config["genome"]["chrom_sizes"],
    output:
        "capcruncher_output/resources/viewpoints/viewpoints.bigBed",
    shell:
        """
        cat {input.viewpoints} | sort -k1,1 -k2,2n > {output}.tmp
        bedToBigBed {output}.tmp {params.chrom_sizes} {output}
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
        comparison=f"[A-Za-z0-9_\.]+-[A-Za-z0-9_\.]+",
        group=f"[A-Za-z0-9_\.]+",
    params:
        color_by=config["hub"].get("color_by", "sample"),
        genome=config["genome"]["name"],
        custom_genome=config["hub"].get("custom_genome", None),
        genome_twobit=config["genome"].get("twobit", None),
        hub_name=config["hub"].get("name"),
        hub_short_label=config["hub"].get("short_label"),
        hub_long_label=config["hub"].get("long_label"),
        hub_email=config["hub"].get("email"),
        genome_organism=config["genome"].get("organism"),
        genome_default_position=config["genome"].get("genome_default_position"),
    script:
        "../scripts/make_ucsc_hub.py"


rule plot:
    input:
        unpack(
            lambda wc: capcruncher.pipeline.utils.get_files_to_plot(
                wc, DESIGN, ASSAY, SAMPLE_NAMES, SUMMARY_METHODS, COMPARE_SAMPLES
            )
        ),
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
        comparison=f"[A-Za-z0-9_\.]+-[A-Za-z0-9_\.]+",
        group=f"[A-Za-z0-9_\.]+",
    log:
        "logs/plot/{viewpoint}.log",
    threads: 1
    script:
        "../scripts/plot.py"


localrules:
    create_ucsc_hub,
    plot,
