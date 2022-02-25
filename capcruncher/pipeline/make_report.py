from audioop import reverse
import os
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import yaml


def plot_deduplication_stats(deduplication_summary_path: os.PathLike):

    df = (
        pd.read_csv(deduplication_summary_path)
        .sort_values("sample")
        .replace("reads_unique", "Unique Reads")
        .replace("reads_removed", "Duplicated Reads")
    )

    fig = px.bar(
        data_frame=df.query('stat_type != "reads_total"'),
        x="stat",
        y="sample",
        color="stat_type",
        template="plotly_white",
        category_orders={
            "sample": sorted(df["sample"].unique()),
            "stat_type": ["Unique Reads", "Duplicated Reads"],
        },
        color_discrete_sequence=["#F9A65A", "grey"],
    )
    # fig.for_each_trace(lambda t: t.update(name=" ".join(t.name.split("_"))))
    fig.update_layout(legend_title_text="")
    fig.update_yaxes(title="")
    fig.update_xaxes(title="Number of Reads")
    fig.update_traces(marker_line_width=0)

    return fig


def plot_trimming_summary(trimming_summary_path: os.PathLike):

    df = pd.read_csv(trimming_summary_path)
    n_samples = len(df["sample"].unique())

    df_summary = df.query(
        'stat_type == "adapters_removed" or stat_type == "reads_total"'
    ).sort_values(["sample", "read_number"])

    subplot_specs = [[{"type": "pie"} for i in range(2)] for j in range(n_samples)]
    fig = make_subplots(
        rows=n_samples,
        cols=2,
        specs=subplot_specs,
        row_titles=sorted(df_summary["sample"].str.replace("_", " ").unique()),
        column_titles=["Read 1", "Read 2"],
    )

    for ii, (sample, df_sample) in enumerate(df_summary.groupby("sample")):
        for jj in range(0, 2):

            df_read_number = df_sample.query(f"read_number == {jj+1}")

            fig.add_trace(
                go.Pie(
                    labels=df_read_number["stat_type"]
                    .str.replace("_", " ")
                    .str.title(),
                    values=df_read_number["stat"],
                    name=f"{sample} {jj+1}",
                    domain={
                        "row": 1,
                    },
                ),
                row=ii + 1,
                col=jj + 1,
            )

    return fig


def format_run_stats_for_flash_figure(
    run_stats_path: os.PathLike,
) -> pd.DataFrame:

    df = pd.read_csv(run_stats_path)
    df_summary = (
        df.loc[df["stage"].isin(["digestion"])]
        .loc[lambda df: df["stat_type"] == "unfiltered"]
        .assign(
            read_type=lambda df: df["read_type"]
            .replace("flashed", "Combined")
            .replace("pe", "Not Combined")
        )
        .groupby(["sample", "stage", "stat_type", "read_type"])["stat"]
        .mean()
        .reset_index()
        .sort_values("sample")
    )

    return df_summary


def plot_flash_summary(run_stats_path: os.PathLike):

    df = format_run_stats_for_flash_figure(run_stats_path)
    fig = px.bar(
        df,
        x="stat",
        y="read_type",
        color="read_type",
        animation_frame="sample",
        range_x=[0, df["stat"].max()],
        template="plotly_white",
        color_discrete_sequence=["#599AD3", "#9E66AB"],
    )

    fig.update_xaxes(title="Number of Read Pairs")
    fig.update_yaxes(title="")
    fig.update_layout(legend_title_text="")
    fig.update_traces(width=0.5)

    try:
        fig["layout"]["updatemenus"] = None
        # fig["layout"]["updatemenus"][0].update(dict(y=1.2, pad={"b": 10, "t": 0, "l": 0}, x=0))
        fig["layout"]["sliders"][0]["pad"] = {"b": 10, "t": 25}
        fig["layout"]["sliders"][0]["x"] = 0
        fig["layout"]["sliders"][0]["len"] = 1

    except (KeyError, IndexError):  # Might only have one sample
        pass

    return fig


def format_digestion_stats_at_read_level(
    digestion_stats_reads_path: os.PathLike,
):

    df = pd.read_csv(digestion_stats_reads_path)

    df = df.query("read_number != 2").assign(
        read_type=lambda df: df["read_type"]
        .replace("flashed", "Combined")
        .replace("pe", "Non-Combined"),
        stat_type=lambda df: df["stat_type"]
        .replace("unfiltered", "All Read Pairs")
        .replace("filtered", "Read Pairs With Valid Slices"),
        sample=lambda df: df["sample"].str.replace("_", " "),
    )
    return df.sort_values("sample")


def plot_digestion_read_summary(digestion_stats_reads_path):

    df = format_digestion_stats_at_read_level(digestion_stats_reads_path)
    fig = px.bar(
        data_frame=df,
        x="stat",
        y="stat_type",
        color="read_type",
        animation_frame="sample",
        template="plotly_white",
        range_x=[0, df.groupby(["sample", "stat_type"])["stat"].sum().max()],
        category_orders={
            "sample": sorted(df["sample"]),
            "read_type": ["Combined", "Non-Combined"],
            "stat_type": ["All Read Pairs", "Read Pairs With Valid Slices"],
        },
        color_discrete_sequence=["#599AD3", "#9E66AB"],
    )
    fig.update_layout(
        legend_title_text="",
        margin={"b": 10},
    )
    fig.update_yaxes(title="")
    # fig.update_xaxes(matches=None, showticklabels=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
    fig.layout["xaxis"]["title"]["text"] = "Number of Slices"
    fig.update_traces(marker_line_width=0)
    fig.update_traces(width=0.5)

    try:
        fig["layout"]["updatemenus"] = None
        # fig["layout"]["updatemenus"][0].update(dict(y=1.2, pad={"b": 10, "t": 0, "l": 0}, x=0))
        fig["layout"]["sliders"][0]["pad"] = {"b": 10, "t": 25}
        fig["layout"]["sliders"][0]["x"] = 0
        fig["layout"]["sliders"][0]["len"] = 1

    except (KeyError, IndexError):  # Might only have one sample
        pass

    return fig


def plot_digestion_histogram(digestion_stats_histogram_path: os.PathLike):

    df = pd.read_csv(digestion_stats_histogram_path)

    df["filtered"] = df["filtered"].map({0: "Pre-filtering", 1: "Post-filtering"})
    df["read_type"] = df["read_type"].map(
        {"flashed": "Combined Reads", "pe": "Non-Combined Reads"}
    )

    df = df.sort_values(["sample", "n_slices", "filtered", "read_type"])

    fig = px.bar(
        data_frame=df,
        x="n_slices",
        y="count",
        color="filtered",
        pattern_shape="read_type",
        animation_frame="sample",
        animation_group="n_slices",
        barmode="group",
        category_orders={"filtered": ["Pre-filtering", "Post-filtering"]},
        template="plotly_white",
        range_x=[0, df["n_slices"].max()],
        range_y=[0, df["count"].max()],
        color_discrete_sequence=["#599AD3", "#9E66AB"],
    )

    fig.update_xaxes(title="")
    fig.update_yaxes(title="")
    fig.update_layout(legend_title_text="")
    fig.update_xaxes(dtick=1)

    try:
        fig["layout"]["updatemenus"] = None
        # fig["layout"]["updatemenus"][0].update(dict(y=1, pad={"b": 10, "t": 0}, x=0))
        fig["layout"]["sliders"][0]["pad"] = {"r": 10, "b": 5, "t": 10}
        fig["layout"]["sliders"][0]["x"] = 0
        fig["layout"]["sliders"][0]["len"] = 1
    except (KeyError, IndexError):
        pass

    return fig


def format_alignment_filtering_read_stats(filtering_read_stats_path: os.PathLike):
    df = pd.read_csv(filtering_read_stats_path)
    df = (
        df.sort_values("stat", ascending=False)
        .query('stat_type != "not-deduplicated"')
        .replace("duplicate_filtered", "partial_duplicate_removal")
        .replace("deduplicated", "full_PCR_duplicate_removal")
        .assign(
            stat_type=lambda df: df["stat_type"]
            .str.replace("_", " ")
            .str.title()
            .str.replace("Pcr", "PCR"),
            read_type=lambda df: df["read_type"]
            .replace("flashed", "Combined")
            .replace("pe", "Non-Combined"),
            sample=lambda df: df["sample"].str.replace("_", " "),
        )
    )
    df.loc[
        (df["stat_type"] == "Full PCR Duplicate Removal")
        & (df["read_type"] == "Non-Combined"),
        "stat",
    ] = (
        df.loc[
            (df["stat_type"] == "Full PCR Duplicate Removal")
            & (df["read_type"] == "Non-Combined"),
            "stat",
        ]
        // 2
    )
    return df.sort_values(
        ["sample", "read_type", "stat"], ascending=[True, True, False]
    )


def plot_alignment_filtering_read_summary(filtering_read_stats_path: os.PathLike):
    df = format_alignment_filtering_read_stats(filtering_read_stats_path)
    fig = px.bar(
        df,
        x="stat",
        y="stat_type",
        color="read_type",
        barmode="group",
        animation_frame="sample",
        animation_group="stat_type",
        template="plotly_white",
        category_orders={
            "stat_type": df["stat_type"].unique(),
            "read_type": list(reversed(["Combined", "Non-Combined"])),
        },
        range_x=[0, df["stat"].max()],
        color_discrete_sequence=list(reversed(["#599AD3", "#9E66AB"])),
    )

    fig.update_xaxes(title="")
    fig.update_yaxes(title="")
    fig.update_layout(legend_title_text="", legend_traceorder="reversed")

    try:
        fig["layout"]["updatemenus"] = None
        # fig["layout"]["updatemenus"][0].update(dict(y=1.1, pad={"b": 5, "t": 0}, x=0))
        fig["layout"]["sliders"][0]["pad"] = {"r": 10, "b": 5, "t": 10}
        fig["layout"]["sliders"][0]["x"] = 0
        fig["layout"]["sliders"][0]["len"] = 1
    except (KeyError, IndexError):
        pass

    return fig


def plot_reporter_summary(reporter_stats_path: os.PathLike):
    df = pd.read_csv(reporter_stats_path)
    df = df.groupby(["sample", "viewpoint", "cis/trans"]).sum().reset_index()
    df = df.replace("cis", "Cis").replace("trans", "Trans")

    fig = px.bar(
        df.sort_values(["sample", "viewpoint"]),
        x="viewpoint",
        y="count",
        color="cis/trans",
        barmode="group",
        animation_frame="sample",
        range_y=[0, df["count"].max()],
        template="plotly_white",
        color_discrete_sequence=["#9CCB86", "#CF597E"],
    )

    fig.update_xaxes(title="")
    fig.update_yaxes(title="")
    fig.update_layout(legend_title_text="")

    try:
        fig["layout"]["updatemenus"] = None
        # fig["layout"]["updatemenus"][0].update(
        #     dict(y=1.2, pad={"l": 0, "b": 10, "t": 0}, x=0)
        # )
        fig["layout"]["sliders"][0]["pad"] = {"r": 0, "b": 5, "t": 50}
        fig["layout"]["sliders"][0]["x"] = 0
        fig["layout"]["sliders"][0]["len"] = 1
    except (KeyError, IndexError):
        pass

    return fig


def format_run_stats_for_overall_summary(run_stats_path: os.PathLike):

    df = pd.read_csv(run_stats_path)
    df = df.sort_values("stat", ascending=False)

    stat_type_mapping = {
        "reads_total": "Total Reads",
        "reads_unique": "PCR Duplicate Filtered (1st pass)",
        "unfiltered": "Passed Trimming and Combining",
        "filtered": "Passed Minimum Slice Length Filter",
        "mapped": "Mapped to Reference genome",
        "contains_single_capture": "Contains one Viewpoint Slice",
        "contains_capture_and_reporter": "Contains one Viewpoint and at least one Reporter Slice",
        "duplicate_filtered": "PCR Duplicate Filtered (2nd pass, partial)",
        "deduplicated": "PCR Duplicate Filtered (final pass)",
    }

    df = df.assign(
        stat_type=lambda df: df["stat_type"].map(stat_type_mapping),
        read_type=lambda df: df["read_type"].replace("flashed", "Combined").replace("pe", "Non-Combined"),
        sample=lambda df: df["sample"].str.replace("_", " "),
    )

    df = df.sort_values("sample")
    return df


def plot_overall_summary(run_stats_path: os.PathLike):

    df = format_run_stats_for_overall_summary(run_stats_path)
    stat_type_order = (
        df.groupby(["sample", "stat_type", "read_type", "read_number"])["stat"]
        .sum()
        .sort_values(ascending=False)
        .reset_index()["stat_type"]
        .unique()
    )

    fig = px.bar(
        df,
        x="stat",
        y="stat_type",
        color="read_type",
        animation_frame="sample",
        animation_group="stat_type",
        barmode="relative",
        template="plotly_white",
        category_orders={
            "sample": sorted(df["sample"].unique()),
            "read_type": ["Combined", "Non-Combined"],
            "stat_type": stat_type_order,
        },
        color_discrete_sequence=["#599AD3", "#9E66AB"],
    )

    fig.update_xaxes(title="")
    fig.update_yaxes(title="")
    fig.update_layout(legend_title_text="")
    fig.update_traces(marker_line_width=0)

    try:
        fig["layout"]["updatemenus"] = None
        # fig["layout"]["updatemenus"][0].update(
        #     dict(y=1.1, pad={"l": 0, "b": 5, "t": 0}, x=0)
        # )
        fig["layout"]["sliders"][0]["pad"] = {"r": 0, "b": 5, "t": 10}
        fig["layout"]["sliders"][0]["x"] = 0
        fig["layout"]["sliders"][0]["len"] = 1
    except (KeyError, IndexError):
        pass

    return fig


def generate_report(
    pipeline_statistics_path: os.PathLike,
    pipeline_report_path: os.PathLike = "report.html",
):

    # Get paths
    fastq_deduplication_path = os.path.join(
        pipeline_statistics_path, "deduplication/deduplication.summary.csv"
    )
    fastq_trimming_path = os.path.join(
        pipeline_statistics_path, "trimming/trimming.summary.csv"
    )
    fastq_digestion_hist_path = os.path.join(
        pipeline_statistics_path, "digestion/digestion.histogram.csv"
    )
    fastq_digestion_read_path = os.path.join(
        pipeline_statistics_path, "digestion/digestion.reads.csv"
    )
    reporter_read_path = os.path.join(
        pipeline_statistics_path, "reporters/reporters.reads.csv"
    )
    reporter_cis_trans_path = os.path.join(
        pipeline_statistics_path, "reporters/reporters.reporters.csv"
    )
    run_stats_path = os.path.join(pipeline_statistics_path, "run_statistics.csv")

    # # Extract HTML template
    # dir_pipeline = os.path.dirname(os.path.abspath(__file__))
    # path_html_template = os.path.join(dir_pipeline, "report_template.html")

    html_header = """
    <html>
    <head>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {
            margin: 0 100;
            background: whitesmoke;
        }
    </style>
    </head>
    <body>
    <h1>Run statistics</h1>
    <p>This report provides statistics for all major pre-processing and filtering steps performed by the pipeline.
        All charts are interactive so hovering over areas of interest will provide additional information.</p>
    """

    html_footer = """</body>
                     </html>"""

    section_template = """    
    <!-- *** Section SECTION_NUMBER *** --->
    <h2>SECTION_NAME</h2>
    <p>SECTION_DESCRIPTION</p>
    FIGURE_HTML
    """
    # <iframe width="1000" height="1000" frameborder="0" seamless="seamless" scrolling="no" \
    # src="PLOT_URL.embed?width=1000&height=1000"></iframe>

    figures = dict(
        deduplication=plot_deduplication_stats(
            fastq_deduplication_path
        ),  # Deduplication
        flashed=plot_flash_summary(run_stats_path),  # Flashed
        digestion_reads=plot_digestion_read_summary(
            fastq_digestion_read_path
        ),  # Digestion reads
        digestion_hist=plot_digestion_histogram(
            fastq_digestion_hist_path
        ),  # Digestion histogram
        alignment_filtering=plot_alignment_filtering_read_summary(
            reporter_read_path
        ),  # Filtering
        reporters=plot_reporter_summary(reporter_cis_trans_path),  # Reporters
        overall=plot_overall_summary(run_stats_path),  # Overall
    )

    figure_name_to_title_mapping = dict(
        deduplication="FASTQ PCR Duplicate Removal",
        flashed="Read pair combination statistics (FLASh)",
        digestion_reads="Fastq <em>in silico</em> digestion statistics (read pair level)",
        digestion_hist="Fastq <em>in silico</em> digestion statistics (slice level)",
        alignment_filtering="Alignment filtering statistics",
        reporters="Identified reporter statistics",
        overall="Pipeline run statistics",
    )

    report_text_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "report_text.yml"
    )
    with open(report_text_path, "r") as r:
        report_text = yaml.safe_load(r)

    with open(pipeline_report_path, "w") as report:

        report.write(html_header)

        for ii, (fig_name, fig) in enumerate(figures.items()):

            fig_html = fig.to_html(
                full_html=False,
                include_plotlyjs=True if ii == 0 else False,
                auto_play=False,
            )
            fig_title = figure_name_to_title_mapping[fig_name]
            fig_text = report_text[fig_title]

            report.write(
                section_template.replace("SECTION_NUMBER", str(ii))
                .replace("SECTION_NAME", fig_title)
                .replace("SECTION_DESCRIPTION", fig_text)
                .replace("FIGURE_HTML", fig_html)
            )

        report.write(html_footer)
