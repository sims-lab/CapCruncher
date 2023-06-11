
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
            "capcruncher_output/pileups/bigwigs/{sample}/{norm}/{sample}_{viewpoint}.bigWig",
            sample=SAMPLE_NAMES,
            norm=["raw", "norm"],
            viewpoint=VIEWPOINT_NAMES,
        ),
        bigwigs_summary=expand(
            "capcruncher_output/comparisons/bigwigs/{group}/{group}.{method}-summary.{viewpoint}.bigWig",
            group=DESIGN["condition"].unique(),
            method=get_summary_methods(),
            viewpoint=VIEWPOINT_NAMES,
        )
        if AGGREGATE_SAMPLES
        else [],
        bigwigs_comparison=expand(
            "capcruncher_output/comparisons/bigwigs/{comparison}/{comparison}.{method}-subtraction.{viewpoint}.bigWig",
            comparison=[
            f"{a}-{b}"
                for a, b in itertools.permutations(DESIGN["condition"].unique(), 2)
            ],
            method=get_summary_methods(),
            viewpoint=VIEWPOINT_NAMES,
        )
        if COMPARE_SAMPLES
        else [],
        report=rules.make_report.output[0],
    output:
        directory(config["hub"]["dir"]),
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


localrules:
    create_ucsc_hub,
