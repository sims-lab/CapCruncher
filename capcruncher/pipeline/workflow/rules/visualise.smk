def get_summary_methods():
    return [
        m
        for m in re.split(r"[,;\s+]", config["compare"].get("summary_methods", "mean,"))
        if m
    ]


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
            "capcruncher_output/comparisons/bigwigs/{group}.{method}-summary.{viewpoint}.bigWig",
            group=DESIGN["condition"].unique(),
            method=get_summary_methods(),
            viewpoint=VIEWPOINT_NAMES,
        )
        if AGGREGATE_SAMPLES
        else [],
        bigwigs_comparison=expand(
            "capcruncher_output/comparisons/bigwigs/{comparison}.{method}-subtraction.{viewpoint}.bigWig",
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


def get_files_to_plot(wc):
    files = {
        "bigwigs": [],
        "subtractions": [],
        "bigwigs_collection": [],
        "heatmaps": [],
    }

    if ASSAY == "tiled":
        files["heatmaps"].extend(
            expand(
                "capcruncher_output/{sample}/{sample}.hdf5",
                sample=SAMPLE_NAMES,
            )
        )
        return files

    if COMPARE_SAMPLES:
        bigwigs_comparison = expand(
            "capcruncher_output/comparisons/bigwigs/{comparison}.{method}-subtraction.{{viewpoint}}.bigWig",
            comparison=[
                f"{a}-{b}"
                for a, b in itertools.permutations(DESIGN["condition"].unique(), 2)
            ],
            method=get_summary_methods(),
        )

        files["subtractions"].extend(bigwigs_comparison)

    bigwigs = expand(
        "capcruncher_output/pileups/bigwigs/{sample}/norm/{sample}_{{viewpoint}}.bigWig",
        sample=SAMPLE_NAMES,
    )

    # if AGGREGATE_SAMPLES:
    #     files["bigwigs_collection"].extend(bigwigs)
    # else:
    #     files["bigwigs"].extend(bigwigs)

    files["bigwigs"].extend(bigwigs)

    return files


def get_plotting_coordinates(wc):
    plot_coords = config["plot"].get("coordinates", None)

    if plot_coords and pathlib.Path(plot_coords).exists():
        df = pd.read_table(
            plot_coords, names=["chrom", "start", "end", "name"], header=None
        )
        df = df.query("name.str.contains(@wc.viewpoint)").iloc[0]

    else:
        df = pd.read_table(
            VIEWPOINTS, names=["chrom", "start", "end", "name"], header=None
        )

        df = df.query("name == @wc.viewpoint").iloc[0]

    return f"{df.chrom}:{df.start}-{df.end}"


rule plot:
    input:
        unpack(get_files_to_plot),
        viewpoints=config["analysis"]["viewpoints"],
    output:
        template="capcruncher_output/figures/{viewpoint}.toml",
        fig="capcruncher_output/figures/{viewpoint}.pdf",
    params:
        coordinates=lambda wc: get_plotting_coordinates(wc),
        viewpoint="{viewpoint}",
        design=DESIGN,
        genes=config["plot"].get("genes", ""),
        binsize=config["analysis"].get("bin_sizes", [None])[0],
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
