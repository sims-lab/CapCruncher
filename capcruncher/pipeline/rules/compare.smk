import itertools


rule union_bedgraph:
    input:
        expand(
            "capcruncher_pileup/04_bedgraph/{sample}/norm/{sample}.norm.{{viewpoint}}.bedgraph",
            sample=DESIGN["sample"].to_list(),
        ),
    output:
        "capcruncher_compare/01_union_bedgraphs/{viewpoint}.tsv",
    params:
        sample_names=" ".join(DESIGN["sample"].to_list()),
    shell:
        """
        bedtools \
        unionbedg \
        -i {input} \
        -header \
        -names {params.sample_names} \
        > {output}
        """


def get_summary_methods():
    return [
        m
        for m in re.split(r"[,;\s+]", config["compare"].get("summary_methods", "mean,"))
        if m
    ]


def identify_columns_based_on_condition():

    condition_args = []
    for group_name, columns in DESIGN.groupby("condition").groups.items():
        condition_args.append(f"-c {','.join(str(c) for c in columns)}")

    condition_args_str = " ".join(condition_args)
    return condition_args_str


rule compare_interactions:
    input:
        "capcruncher_compare/01_union_bedgraphs/{viewpoint}.tsv",
    output:
        bedgraphs_summary=expand(
            "capcruncher_compare/02_compare_interactions/{group}.{method}-summary.{{viewpoint}}.bedgraph",
            group=DESIGN["condition"].unique(),
            method=get_summary_methods(),
        ),
        bedgraphs_compare=expand(
            "capcruncher_compare/02_compare_interactions/{comparison}.{method}-subtraction.{{viewpoint}}.bedgraph",
            comparison=[
            f"{a}-{b}"
                for a, b in itertools.permutations(DESIGN["condition"].unique(), 2)
            ],
            method=get_summary_methods(),
        ),
    params:
        output_prefix="capcruncher_compare/02_compare_interactions/",
        summary_methods=" ".join([f"-m {m}" for m in get_summary_methods()]),
        names=" ".join([f"-n {group}" for group in DESIGN["condition"].unique()]),
        conditions=identify_columns_based_on_condition(),
    log:
        "logs/compare_interactions/{viewpoint}.log",
    shell:
        """
        capcruncher \
        interactions \
        compare \
        summarise \
        {input} \
        -o {params.output_prefix} \
        -f bedgraph \
        {params.summary_methods} \
        {params.names} \
        {params.conditions} \
        --subtraction \
        --suffix .{wildcards.viewpoint} \
        > {log} 2>&1
        """


use rule bedgraph_to_bigwig as bigwig_compared with:
    input:
        bedgraph="capcruncher_compare/02_compare_interactions/{comparison}.{method}-subtraction.{viewpoint}.bedgraph",
    output:
        bigwig="capcruncher_compare/03_bigwig/{comparison}/{comparison}.{method}-subtraction.{viewpoint}.bw",
    params:
        chrom_sizes=config["genome"]["chrom_sizes"],
    log:
        "logs/bedgraph_to_bigwig/{comparison}.{method}-subtraction.{viewpoint}.log",


use rule bedgraph_to_bigwig as bigwig_summarised with:
    input:
        bedgraph="capcruncher_compare/02_compare_interactions/{group}.{method}-summary.{viewpoint}.bedgraph",
    output:
        bigwig="capcruncher_compare/03_bigwig/{group}/{group}.{method}-summary.{viewpoint}.bw",
    params:
        chrom_sizes=config["genome"]["chrom_sizes"],
    log:
        "logs/bedgraph_to_bigwig/{group}.{method}-summary.{viewpoint}.log",
