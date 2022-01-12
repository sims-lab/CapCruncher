import os
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def plot_deduplication_stats(deduplication_summary_path: os.PathLike):

    df = pd.read_csv(deduplication_summary_path)

    fig = px.bar(
        data_frame=df.query('stat_type != "reads_total"'),
        x="stat",
        y="sample",
        color="stat_type",
        template="simple_white",
        category_orders={
            "sample": sorted(df["sample"].unique()),
            "stat_type": ("reads_unique", "reads_removed"),
        },
        color_discrete_sequence=["#1f77b4", "grey"],
    )
    fig.for_each_trace(lambda t: t.update(name=" ".join(t.name.split("_"))))
    fig.update_layout(legend_title_text="")
    fig.update_yaxes(title="Sample")
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
            .replace("flashed", "Flashed")
            .replace("pe", "PE")
        )
        .groupby(["sample", "stage", "stat_type", "read_type"])["stat"]
        .mean()
        .reset_index()
    )

    return df_summary


def plot_flash_summary(run_stats_path: os.PathLike):

    df = format_run_stats_for_flash_figure(run_stats_path)
    fig = px.bar(
        data_frame=df,
        x="stat",
        y="sample",
        color="read_type",
        template="simple_white",
        category_orders={
            "sample": sorted(df["sample"]),
            "read_type": ["Flashed", "PE"],
        },
    )
    fig.update_layout(legend_title_text="")
    fig.update_yaxes(title="Sample")
    fig.update_xaxes(title="Number of Read Pairs")
    fig.update_traces(marker_line_width=0)

    return fig


def format_digestion_stats_at_read_level(
    digestion_stats_reads_path: os.PathLike,
):

    df = pd.read_csv(digestion_stats_reads_path)

    df = df.query("read_number != 2").assign(
        read_type=lambda df: df["read_type"]
        .replace("flashed", "Flashed")
        .replace("pe", "PE"),
        stat_type=lambda df: df["stat_type"]
        .replace("unfiltered", "All Read Pairs")
        .replace("filtered", "Reads with slices"),
        sample=lambda df: df["sample"].str.replace("_", " "),
    )
    return df


def plot_digestion_read_summary(digestion_stats_reads_path):

    df = format_digestion_stats_at_read_level(digestion_stats_reads_path)
    fig = px.bar(
        data_frame=df,
        x="stat",
        y="stat_type",
        color="read_type",
        facet_row="sample",
        template="simple_white",
        category_orders={
            "sample": sorted(df["sample"]),
            "read_type": ["Flashed", "PE"],
        },
    )
    fig.update_layout(
        legend_title_text="",
        margin={"b": 10},
    )
    fig.update_yaxes(title="", autorange="reversed")
    fig.update_xaxes(matches=None, showticklabels=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
    fig.layout["xaxis"]["title"]["text"] = "Number of Slices (Reads with RE sites)"
    fig.update_traces(marker_line_width=0)

    return fig


def plot_digestion_histogram(digestion_stats_histogram_path: os.PathLike):

    df = pd.read_csv(digestion_stats_histogram_path)

    fig = px.histogram(
        data_frame=df.assign(
            read_number=lambda df: df["read_number"].map(
                {0: "Flashed", 1: "PE R1", 2: "PE R2"}
            )
        ),
        x="count",
        y="n_slices",
        color="read_number",
        facet_row="sample",
        template="simple_white",
        barmode="group",
        hover_data=["count"],
        category_orders={"read_number": ["Flashed", "PE R1", "PE R2"]},
    )

    fig.update_layout(legend_title_text="")
    fig.update_yaxes(title="Frequency", matches=None, showticklabels=True)
    fig.update_xaxes(dtick=1, showticklabels=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
    fig.update_traces(marker_line_width=0)

    return fig


def format_alignment_filtering_read_stats(filtering_read_stats_path: os.PathLike):
    df = pd.read_csv(filtering_read_stats_path)
    return (
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
            .replace("flashed", "Flashed")
            .replace("pe", "PE"),
            sample=lambda df: df["sample"].str.replace("_", " "),
        )
    )


def plot_alignment_filtering_read_summary(filtering_read_stats_path: os.PathLike):
    df = format_alignment_filtering_read_stats(filtering_read_stats_path)
    fig = px.bar(
        data_frame=df.sort_values("stat", ascending=False),
        x="stat",
        y="stat_type",
        template="simple_white",
        color="read_type",
        facet_row="sample",
        category_orders={
            "stat_type": df["stat_type"].unique(),
            "read_type": ["Flashed", "PE"],
            "sample": sorted(df["sample"].unique()),
        },
    )
    fig.update_xaxes(matches=None, showticklabels=True)
    fig.update_yaxes(title="")
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
    fig.update_layout(legend_title_text="")
    fig.update_traces(marker_line_width=0)

    return fig


def plot_reporter_summary(reporter_stats_path: os.PathLike):
    df = pd.read_csv(reporter_stats_path)
    n_probes = df["viewpoint"].nunique()

    fig = px.bar(
        data_frame=df.groupby(["sample", "viewpoint", "cis/trans"])
        .agg({"count": "sum"})
        .reset_index()
        .assign(sample=lambda df: df["sample"].str.replace("_", " ")),
        x="count",
        y="viewpoint",
        color="cis/trans",
        facet_row="sample",
        barmode="group",
        template="simple_white",
        category_orders={
            "cis/trans": ["trans", "cis"],
            "viewpoint": sorted(df["viewpoint"].unique()),
            "sample": sorted(df["sample"].unique()),
        },
        labels={"count": "Number of reporters"},
    )
    fig.update_yaxes(title_text="")
    fig.update_xaxes(matches=None, showticklabels=True)
    fig.for_each_trace(lambda t: t.update(name=t.name.split("=")[0]))
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
    fig.update_layout(legend={"traceorder": "reversed", "title": ""})
    fig.update_traces(marker_line_width=0)
    return fig
