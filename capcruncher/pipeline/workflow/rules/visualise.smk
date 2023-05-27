def get_hub_pileup_files(wildcards):

    bigwigs = {}
    bigwigs["bigwigs"] = expand(
        "capcruncher_output/pileups/bigwigs/{sample}/{norm}/{sample}_{viewpoint}.bigWig",
        sample=SAMPLE_NAMES,
        norm=["raw", "norm"],
        viewpoint=get_existing_viewpoints(wildcards),
    )

    if COMPARE_SAMPLES:

        bigwigs["bigwigs_compared"] = expand(
            "capcruncher_output/comparisons/bigwigs/{comparison}/{comparison}.{method}-subtraction.{viewpoint}.bigWig",
            comparison=[
                f"{a}-{b}"
                for a, b in itertools.permutations(DESIGN["condition"].unique(), 2)
            ],
            method=get_summary_methods(),
            viewpoint=get_existing_viewpoints,
        )
        bigwigs["bigwigs_summarised"] = expand(
            "capcruncher_output/comparisons/bigwigs/{group}/{group}.{method}-summary.{viewpoint}.bigWig",
            group=DESIGN["condition"].unique(),
            method=get_summary_methods(),
            viewpoint=get_existing_viewpoints(wildcards),
        )

    return bigwigs


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
        bigwigs=get_hub_pileup_files,
        viewpoints=rules.viewpoints_to_bigbed.output[0],
        report=rules.make_report.output[0],
    output:
        directory(config["hub"]["dir"]),
    script:
        "../scripts/make_ucsc_hub.py"


localrules:
    create_ucsc_hub,
