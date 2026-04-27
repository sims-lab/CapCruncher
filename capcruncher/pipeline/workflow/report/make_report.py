import html
import json
import pathlib
from collections.abc import Iterable

import pandas as pd
import plotly.express as px
import yaml


REPORT_TITLE = "CapCruncher Run Report"
COLORWAY = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"]
READ_TYPE_COLORS = ["#0072B2", "#E69F00"]
DEDUP_COLORS = ["#009E73", "#6C757D"]
COMMON_LABELS = {
    "count": "Count",
    "filter_status": "Digestion Status",
    "n_fragments": "Fragments",
    "n_reads": "Reads",
    "n_slices": "Slices",
    "read_number": "Read",
    "read_type": "Read Type",
    "sample": "Sample",
    "stage": "Stage",
    "stage_label": "Stage",
    "viewpoint": "Viewpoint",
    "cis_or_trans": "Reporter Type",
}
STAGE_LABELS = {
    "pre-filtering": "Pre-filtering",
    "mapped": "Mapped",
    "contains_single_capture": "Contains One Viewpoint",
    "contains_capture_and_reporter": "Contains Viewpoint and Reporter",
    "duplicate_filtered": "Partial PCR Duplicate Removal",
    "final_duplicate_removal": "Final PCR Duplicate Removal",
}


def load_json(path: str | pathlib.Path):
    with pathlib.Path(path).open() as handle:
        return json.load(handle)


def load_json_strings(path: str | pathlib.Path) -> list[dict]:
    entries = load_json(path)
    return [json.loads(entry) if isinstance(entry, str) else entry for entry in entries]


def natural_sort_paths(paths: Iterable[str | pathlib.Path]) -> list[pathlib.Path]:
    return sorted(pathlib.Path(path) for path in paths)


def require_paths(paths: Iterable[str | pathlib.Path], label: str) -> list[pathlib.Path]:
    path_list = natural_sort_paths(paths)
    if not path_list:
        raise ValueError(f"No {label} statistics files were found for the report")
    return path_list


def read_type_label(value: str) -> str:
    labels = {
        "flashed": "Combined",
        "Flashed": "Combined",
        "pe": "Non-Combined",
        "Pe": "Non-Combined",
        "unflashed": "Non-Combined",
    }
    return labels.get(value, value)


def read_number_label(value: str) -> str:
    labels = {"read1": "Read 1", "read2": "Read 2"}
    return labels.get(value, value)


def stage_label(value: str) -> str:
    return STAGE_LABELS.get(value, value.replace("_", " ").title())


def reporter_type_label(value: str) -> str:
    labels = {"cis": "Cis", "trans": "Trans"}
    return labels.get(value, value)


def frame_to_table_html(df: pd.DataFrame) -> str:
    table = df.to_html(
        classes="data-table",
        index=False,
        border=0,
        escape=False,
        table_id=None,
        float_format=lambda value: f"{value:.2f}",
    )
    return table.replace("<table ", '<table data-interactive="true" ')


def plot_html(fig, include_plotlyjs: bool) -> str:
    return fig.to_html(
        full_html=False,
        include_plotlyjs=True if include_plotlyjs else False,
        auto_play=False,
        config={"responsive": True, "displaylogo": False},
    )


def polish_figure(
    fig,
    x_title: str | None = None,
    y_title: str | None = None,
    x_tickangle: int | None = None,
):
    fig.update_layout(
        colorway=COLORWAY,
        font={"family": "-apple-system, BlinkMacSystemFont, Segoe UI, sans-serif"},
        hoverlabel={
            "bgcolor": "#ffffff",
            "bordercolor": "#d7d9d2",
            "font_size": 13,
        },
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "left",
            "x": 0,
            "title_text": "",
        },
        margin={"l": 150, "r": 36, "t": 78, "b": 96},
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
    )
    fig.update_xaxes(
        title=x_title,
        automargin=True,
        gridcolor="#edf0eb",
        linecolor="#d7d9d2",
        tickangle=x_tickangle,
        ticklabeloverflow="allow",
        zerolinecolor="#d7d9d2",
    )
    fig.update_yaxes(
        title=y_title,
        automargin=True,
        gridcolor="#edf0eb",
        linecolor="#d7d9d2",
        ticklabeloverflow="allow",
        zerolinecolor="#d7d9d2",
    )
    fig.update_traces(marker_line_width=0, selector={"type": "bar"})
    fig.update_traces(line_width=2.5, selector={"type": "scatter"})
    return fig


def tab_group(*items: tuple[str, str]) -> str:
    buttons = []
    panels = []
    for index, (title, body) in enumerate(items):
        selected = "true" if index == 0 else "false"
        hidden = "" if index == 0 else " hidden"
        buttons.append(
            f'<button type="button" role="tab" aria-selected="{selected}">'
            f"{html.escape(title)}</button>"
        )
        panels.append(f'<div role="tabpanel"{hidden}>{body}</div>')

    return (
        '<div class="tabs">'
        f'<div class="tab-list" role="tablist">{"".join(buttons)}</div>'
        f'{"".join(panels)}'
        "</div>"
    )


def figure_range_max(series: pd.Series) -> int | float:
    maximum = series.max()
    if pd.isna(maximum) or maximum == 0:
        return 1
    return maximum + maximum * 0.1


def load_fastq_deduplication(paths: Iterable[str | pathlib.Path]) -> pd.DataFrame:
    df = pd.DataFrame([load_json(path) for path in paths])
    return df.assign(
        unique=lambda frame: frame["total"] - frame["duplicates"],
        percentage=lambda frame: frame["duplicates"] / frame["total"] * 100,
    )


def load_fastq_trimming(path: str | pathlib.Path) -> pd.DataFrame:
    df = pd.DataFrame(load_json_strings(path))
    return df.assign(
        percentage_trimmed=lambda frame: frame.get(
            "percentage_trimmed",
            frame["reads_with_adapter_identified"] / frame["reads_input"] * 100,
        ),
        percentage_passing_quality_filter=lambda frame: frame.get(
            "percentage_passing_quality_filter",
            frame["reads_output"] / frame["reads_input"] * 100,
        ),
    )


def load_fastq_flash(path: str | pathlib.Path) -> pd.DataFrame:
    df = pd.DataFrame(load_json_strings(path))
    return df.assign(
        n_total=lambda frame: frame.get(
            "n_total", frame["n_combined"] + frame["n_uncombined"]
        ),
        percentage_combined=lambda frame: frame.get(
            "percentage_combined",
            frame["n_combined"] / (frame["n_combined"] + frame["n_uncombined"]) * 100,
        ),
    )


def load_digestion(paths: Iterable[str | pathlib.Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    histograms = []

    for stat in [load_json(path) for path in paths]:
        for filter_status in ["unfiltered", "filtered"]:
            read_pair_stat = stat["read_stats"][filter_status]
            for read_number in ["read1", "read2"]:
                count = read_pair_stat.get(read_number)
                if count is None:
                    continue
                rows.append(
                    {
                        "sample": stat["sample"],
                        "read_type": read_type_label(stat["read_type"]),
                        "filter_status": filter_status.replace(
                            "unfiltered", "Pre-digestion"
                        ).replace("filtered", "Post-digestion"),
                        "read_number": read_number_label(read_number),
                        "count": count,
                    }
                )

        length_frame = histogram_pair_to_dataframe(stat["histograms"]["lengths"])
        if not length_frame.empty:
            histograms.append(
                length_frame.assign(
                    sample=stat["sample"],
                    read_type=read_type_label(stat["read_type"]),
                )
            )

    digestion_reads = pd.DataFrame(rows)
    length_histograms = pd.concat(histograms, ignore_index=True)
    return digestion_reads, length_histograms


def histogram_pair_to_dataframe(read_pair: dict) -> pd.DataFrame:
    frames = []
    for read_number in ["read1", "read2"]:
        histogram = read_pair.get(read_number)
        if not histogram:
            continue
        name = histogram.get("name", "value")
        frame = pd.DataFrame(
            [(int(value), count) for value, count in histogram.get("hist", {}).items()],
            columns=[name, "count"],
        )
        if not frame.empty:
            frames.append(frame.assign(read_number=read_number_label(read_number)))
    if not frames:
        return pd.DataFrame(columns=["slice_length", "count", "read_number"])
    return pd.concat(frames, ignore_index=True)


def load_filtering(
    filtering_paths: Iterable[str | pathlib.Path],
    deduplication_paths: Iterable[str | pathlib.Path],
) -> pd.DataFrame:
    filter_stats = [load_json(path) for path in filtering_paths]
    rows = [stat for stats in filter_stats for stat in stats["stats"]]
    df_filter_stats = pd.DataFrame(rows)

    stage_order = {
        stage: order
        for order, stage in enumerate(df_filter_stats["stage"].drop_duplicates())
    }
    df_filter_stats = (
        df_filter_stats.groupby(["sample", "stage", "read_type"], as_index=False)
        .sum()
        .assign(
            read_type=lambda df: df["read_type"].map(read_type_label),
            stage_order=lambda df: df["stage"].map(stage_order),
        )
    )

    dedup_stats = [load_json(path) for path in deduplication_paths]
    if dedup_stats:
        df_dedup = pd.DataFrame(dedup_stats)
        df_filter_stats = pd.concat(
            [
                df_filter_stats,
                df_dedup[
                    ["sample", "read_type", "n_unique_reads", "n_unique_slices"]
                ]
                .rename(
                    columns={
                        "n_unique_reads": "n_fragments",
                        "n_unique_slices": "n_slices",
                    }
                )
                .assign(
                    read_type=lambda df: df["read_type"].map(read_type_label),
                    stage="final_duplicate_removal",
                    stage_order=df_filter_stats["stage_order"].max() + 1,
                ),
            ],
            ignore_index=True,
        )

    return df_filter_stats.sort_values(["sample", "stage_order", "read_type"])


def load_cis_trans(paths: Iterable[str | pathlib.Path]) -> pd.DataFrame:
    stats = [load_json(path) for path in paths]
    rows = [stat for sample in stats for stat in sample["stats"]]
    return (
        pd.DataFrame(rows)
        .assign(
            read_type=lambda df: df["read_type"].map(read_type_label),
            cis_or_trans=lambda df: df["cis_or_trans"].map(reporter_type_label),
        )
        .groupby(["sample", "read_type", "cis_or_trans", "viewpoint"], as_index=False)
        .sum()
    )


def build_report(
    output: str | pathlib.Path,
    fastq_deduplication_paths: Iterable[str | pathlib.Path],
    fastq_trimming_path: str | pathlib.Path,
    fastq_flash_path: str | pathlib.Path,
    fastq_digestion_paths: Iterable[str | pathlib.Path],
    reporter_filtering_paths: Iterable[str | pathlib.Path],
    reporter_deduplication_paths: Iterable[str | pathlib.Path],
    reporter_cis_trans_paths: Iterable[str | pathlib.Path],
):
    df_dedup = load_fastq_deduplication(
        require_paths(fastq_deduplication_paths, "FASTQ deduplication")
    )
    df_trim = load_fastq_trimming(fastq_trimming_path)
    df_flash = load_fastq_flash(fastq_flash_path)
    df_digestion, df_lengths = load_digestion(
        require_paths(fastq_digestion_paths, "FASTQ digestion")
    )
    df_filter = load_filtering(
        require_paths(reporter_filtering_paths, "reporter filtering"),
        require_paths(reporter_deduplication_paths, "reporter deduplication"),
    )
    df_cis_trans = load_cis_trans(
        require_paths(reporter_cis_trans_paths, "cis/trans reporter")
    )

    sections = make_sections(
        df_dedup=df_dedup,
        df_trim=df_trim,
        df_flash=df_flash,
        df_digestion=df_digestion,
        df_lengths=df_lengths,
        df_filter=df_filter,
        df_cis_trans=df_cis_trans,
    )

    pathlib.Path(output).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(output).write_text(render_html(sections), encoding="utf-8")


def make_sections(
    df_dedup: pd.DataFrame,
    df_trim: pd.DataFrame,
    df_flash: pd.DataFrame,
    df_digestion: pd.DataFrame,
    df_lengths: pd.DataFrame,
    df_filter: pd.DataFrame,
    df_cis_trans: pd.DataFrame,
) -> list[dict[str, str]]:
    df_dedup_summary = (
        df_dedup.groupby("sample", as_index=False).sum(numeric_only=True).assign(
            percentage=lambda df: df["duplicates"] / df["total"] * 100,
            unique=lambda df: df["total"] - df["duplicates"],
        )
    )

    sections = []
    include_plotlyjs = True

    df = df_dedup_summary[["sample", "duplicates", "unique"]].melt(
        id_vars="sample", var_name="read_type", value_name="count"
    ).assign(
        read_type=lambda frame: frame["read_type"].map(
            {"unique": "Unique Reads", "duplicates": "Duplicate Reads"}
        )
    )
    fig = px.bar(
        df,
        x="count",
        y="sample",
        color="read_type",
        template="plotly_white",
        labels=COMMON_LABELS,
        category_orders={"read_type": ["Unique Reads", "Duplicate Reads"]},
        color_discrete_sequence=DEDUP_COLORS,
    )
    polish_figure(fig, x_title="Reads", y_title="")
    sections.append(
        section(
            "FASTQ PCR Duplicate Removal",
            tab_group(
                (
                    "Table",
                    frame_to_table_html(
                        df_dedup_summary[
                            ["sample", "total", "unique", "duplicates", "percentage"]
                        ]
                        .round(2)
                        .rename(
                            columns={
                                "sample": "Sample",
                                "total": "Total Reads",
                                "unique": "Unique Reads",
                                "duplicates": "Duplicate Reads",
                                "percentage": "Percentage Duplicated",
                            }
                        )
                    ),
                ),
                ("Bar Chart", plot_html(fig, include_plotlyjs)),
            ),
        )
    )
    include_plotlyjs = False

    sections.append(
        section(
            "Trimming",
            frame_to_table_html(
                df_trim.rename(columns=lambda col: col.replace("_", " ").title()).round(
                    2
                )
            ),
        )
    )

    sections.append(
        section(
            "Read pair combination statistics (FLASh)",
            frame_to_table_html(
                df_flash.rename(columns=lambda col: col.replace("_", " ").title()).round(
                    2
                )
            ),
        )
    )

    df_digestion_table = (
        df_digestion.pivot_table(
            index=["sample", "read_type", "read_number"],
            columns="filter_status",
            values="count",
            aggfunc="sum",
            fill_value=0,
        )
        .reset_index()
        .assign(
            percentage_digested=lambda df: (
                df["Post-digestion"] / df["Pre-digestion"] * 100
            )
            .replace([float("inf"), -float("inf")], 0)
            .fillna(0)
        )
    )
    df_digestion_table.columns.name = None
    fig = px.line(
        df_digestion,
        x="filter_status",
        y="count",
        color="read_number",
        line_dash="read_type",
        animation_frame="sample",
        markers=True,
        range_y=(0, figure_range_max(df_digestion["count"])),
        template="plotly_white",
        labels=COMMON_LABELS,
        category_orders={
            "filter_status": ["Pre-digestion", "Post-digestion"],
            "read_number": ["Read 1", "Read 2"],
            "read_type": ["Combined", "Non-Combined"],
        },
        color_discrete_sequence=READ_TYPE_COLORS,
    )
    polish_figure(fig, x_title="Digestion Status", y_title="Reads")
    strip_animation_controls(fig)
    sections.append(
        section(
            "Fastq <em>in silico</em> digestion statistics (read pair level)",
            tab_group(
                (
                    "Table",
                    frame_to_table_html(
                        df_digestion_table.rename(
                            columns=lambda col: col.title().replace("_", " ")
                        ).round(2)
                    ),
                ),
                ("Digestion Filtering Plot", plot_html(fig, include_plotlyjs)),
            ),
        )
    )

    fig = px.histogram(
        df_lengths.sort_values("sample"),
        x="slice_length",
        y="count",
        pattern_shape="read_number",
        facet_col="read_type",
        color="read_type",
        animation_frame="sample",
        range_x=(0, df_lengths["slice_length"].max() + 1),
        barmode="group",
        template="plotly_white",
        labels=COMMON_LABELS | {"slice_length": "Slice Length"},
        category_orders={
            "read_number": ["Read 1", "Read 2"],
            "read_type": ["Combined", "Non-Combined"],
        },
        color_discrete_sequence=READ_TYPE_COLORS,
    )
    fig.for_each_annotation(
        lambda annotation: annotation.update(text=annotation.text.split("=")[-1])
    )
    polish_figure(fig, x_title="Slice Length", y_title="Slices")
    strip_animation_controls(fig)
    sections.append(
        section(
            "Fastq <em>in silico</em> digestion statistics (slice level)",
            plot_html(fig, include_plotlyjs),
        )
    )

    df_filter_plot = df_filter.assign(
        stage_label=lambda df: df["stage"].map(stage_label)
    )
    fig_fragments = px.line(
        df_filter_plot,
        x="stage_label",
        y="n_fragments",
        color="read_type",
        animation_frame="sample",
        markers=True,
        range_y=(0, figure_range_max(df_filter_plot["n_fragments"])),
        template="plotly_white",
        labels=COMMON_LABELS,
        category_orders={"read_type": ["Combined", "Non-Combined"]},
        color_discrete_sequence=READ_TYPE_COLORS,
    )
    polish_figure(
        fig_fragments,
        x_title="Filtering Stage",
        y_title="Fragments",
        x_tickangle=-30,
    )
    strip_animation_controls(fig_fragments)

    fig_slices = px.line(
        df_filter_plot,
        x="stage_label",
        y="n_slices",
        color="read_type",
        animation_frame="sample",
        markers=True,
        range_y=(0, figure_range_max(df_filter_plot["n_slices"])),
        template="plotly_white",
        labels=COMMON_LABELS,
        category_orders={"read_type": ["Combined", "Non-Combined"]},
        color_discrete_sequence=READ_TYPE_COLORS,
    )
    polish_figure(
        fig_slices,
        x_title="Filtering Stage",
        y_title="Slices",
        x_tickangle=-30,
    )
    strip_animation_controls(fig_slices)
    sections.append(
        section(
            "Alignment filtering statistics",
            tab_group(
                (
                    "Table",
                    frame_to_table_html(
                        df_filter_plot[
                            [
                                "sample",
                                "read_type",
                                "stage_label",
                                "n_fragments",
                                "n_slices",
                            ]
                        ].rename(
                            columns={
                                "sample": "Sample",
                                "read_type": "Read Type",
                                "stage_label": "Stage",
                                "n_fragments": "Fragments",
                                "n_slices": "Slices",
                            }
                        )
                    ),
                ),
                ("Fragments", plot_html(fig_fragments, include_plotlyjs)),
                ("Slices", plot_html(fig_slices, include_plotlyjs)),
            ),
        )
    )

    fig = px.bar(
        df_cis_trans,
        x="count",
        y="viewpoint",
        color="cis_or_trans",
        pattern_shape="read_type",
        animation_frame="sample",
        facet_row="cis_or_trans",
        template="plotly_white",
        labels=COMMON_LABELS,
        category_orders={
            "cis_or_trans": ["Cis", "Trans"],
            "read_type": ["Combined", "Non-Combined"],
        },
        color_discrete_sequence=READ_TYPE_COLORS,
        range_x=(0, figure_range_max(df_cis_trans["count"])),
    )
    fig.for_each_annotation(
        lambda annotation: annotation.update(text=annotation.text.split("=")[-1])
    )
    polish_figure(fig, x_title="Reporters", y_title="Viewpoint")
    strip_animation_controls(fig)
    sections.append(
        section(
            "Identified reporter statistics",
            frame_to_table_html(
                df_cis_trans.rename(
                    columns={
                        "sample": "Sample",
                        "read_type": "Read Type",
                        "cis_or_trans": "Cis/Trans",
                        "viewpoint": "Viewpoint",
                        "count": "Count",
                    }
                )
            )
            + plot_html(fig, include_plotlyjs),
        )
    )

    df_stats = make_overall_stats(df_dedup_summary, df_trim, df_digestion, df_filter)
    fig = px.line(
        df_stats,
        x="stage",
        y="n_reads",
        animation_frame="sample",
        template="plotly_white",
        range_y=(0, figure_range_max(df_stats["n_reads"])),
        markers=True,
        labels=COMMON_LABELS,
    )
    polish_figure(fig, x_title="", y_title="Reads", x_tickangle=-30)
    strip_animation_controls(fig)
    sections.append(
        section(
            "Pipeline run statistics",
            frame_to_table_html(
                df_stats[["sample", "stage", "n_reads"]].rename(
                    columns={
                        "sample": "Sample",
                        "stage": "Stage",
                        "n_reads": "Reads",
                    }
                )
            )
            + plot_html(fig, include_plotlyjs),
        )
    )

    return sections


def make_overall_stats(
    df_dedup: pd.DataFrame,
    df_trim: pd.DataFrame,
    df_digestion: pd.DataFrame,
    df_filter: pd.DataFrame,
) -> pd.DataFrame:
    raw = df_dedup[["sample", "total"]].rename(columns={"total": "n_reads"}).assign(
        stage="Raw", stage_order=0
    )
    fastq_dedup = df_dedup[["sample", "unique"]].rename(
        columns={"unique": "n_reads"}
    ).assign(stage="Fastq Deduplication", stage_order=1)
    fastq_trim = (
        df_trim[["sample", "reads_output"]]
        .drop_duplicates()
        .groupby("sample", as_index=False)
        .min()
        .rename(columns={"reads_output": "n_reads"})
        .assign(stage="Fastq Trimming", stage_order=2)
    )
    fastq_digest = (
        df_digestion.query(
            "filter_status == 'Post-digestion' and read_number == 'Read 1'"
        )
        .groupby("sample", as_index=False)["count"]
        .sum()
        .rename(columns={"count": "n_reads"})
        .assign(stage="Fastq Digestion", stage_order=3)
    )
    aln_filter = (
        df_filter.groupby(["sample", "stage", "stage_order"], as_index=False)[
            "n_fragments"
        ]
        .sum()
        .rename(columns={"n_fragments": "n_reads"})
        .assign(
            stage=lambda df: df["stage"].str.replace("_", " ").str.title(),
            stage_order=lambda df: df["stage_order"] + 4,
        )
    )
    return pd.concat(
        [raw, fastq_dedup, fastq_trim, fastq_digest, aln_filter], ignore_index=True
    ).sort_values(["sample", "stage_order"])


def strip_animation_controls(fig):
    fig.update_layout(updatemenus=[])
    if fig.layout.sliders:
        fig.layout.sliders[0].currentvalue = {"prefix": "Sample: "}
        fig.layout.sliders[0].pad = {"r": 10, "b": 5, "t": 10}
        fig.layout.sliders[0].x = 0
        fig.layout.sliders[0].len = 1


def section(title: str, body: str) -> dict[str, str]:
    return {"title": title, "body": body}


def render_html(sections: list[dict[str, str]]) -> str:
    text_path = pathlib.Path(__file__).with_name("report_text.yml")
    report_text = yaml.safe_load(text_path.read_text(encoding="utf-8"))

    nav = "\n".join(
        f'<a href="#section-{idx}">{html.escape(strip_tags(item["title"]))}</a>'
        for idx, item in enumerate(sections)
    )
    rendered_sections = "\n".join(
        render_section(idx, item, report_text) for idx, item in enumerate(sections)
    )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{REPORT_TITLE}</title>
<style>
html {{
  scroll-behavior: smooth;
}}
body {{
  margin: 0;
  color: #23272a;
  background: #f4f5f2;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  line-height: 1.45;
}}
header {{
  background: #ffffff;
  border-bottom: 1px solid #d7d9d2;
  padding: 30px 40px 22px;
}}
main {{
  max-width: 1220px;
  margin: 0 auto;
  padding: 30px 24px 56px;
}}
h1 {{
  margin: 0 0 8px;
  font-size: 34px;
  letter-spacing: 0;
}}
h2 {{
  margin-top: 0;
  font-size: 23px;
  letter-spacing: 0;
}}
nav {{
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  margin-top: 20px;
}}
nav a {{
  border: 1px solid #c8d3d8;
  border-radius: 999px;
  color: #245a70;
  background: #f8fbfb;
  padding: 5px 10px;
  text-decoration: none;
  font-size: 14px;
}}
nav a:hover {{
  border-color: #245a70;
  background: #eef6f6;
}}
section {{
  background: #fff;
  border: 1px solid #dfe2dc;
  border-radius: 8px;
  box-shadow: 0 1px 2px rgba(22, 28, 31, 0.05);
  margin: 0 0 26px;
  padding: 26px;
}}
.description {{
  max-width: 900px;
  color: #4d5356;
}}
.table-wrap {{
  max-width: 100%;
  overflow-x: auto;
  margin: 18px 0;
}}
.table-filter {{
  width: min(360px, 100%);
  box-sizing: border-box;
  border: 1px solid #c8cfca;
  border-radius: 6px;
  padding: 9px 11px;
  margin: 12px 0 10px;
  font-size: 14px;
  background: #fff;
}}
.data-table {{
  border-collapse: collapse;
  width: 100%;
  font-size: 14px;
}}
.data-table th,
.data-table td {{
  border-bottom: 1px solid #e5e8e2;
  padding: 9px 11px;
  text-align: left;
  white-space: nowrap;
}}
.data-table th {{
  background: #eef3f1;
  font-weight: 650;
  cursor: pointer;
  user-select: none;
}}
.data-table tbody tr:nth-child(even) {{
  background: #fafbf9;
}}
.data-table tbody tr:hover {{
  background: #f2f7f6;
}}
.plotly-graph-div {{
  margin-top: 14px;
}}
.tab-list {{
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  border-bottom: 1px solid #d7d9d2;
  margin: 18px 0;
}}
.tab-list button {{
  appearance: none;
  border: 1px solid #d7d9d2;
  border-bottom: 0;
  border-radius: 6px 6px 0 0;
  background: #f4f5f2;
  color: #293033;
  padding: 9px 13px;
  font: inherit;
  cursor: pointer;
}}
.tab-list button[aria-selected="true"] {{
  background: #fff;
  color: #245a70;
  font-weight: 650;
}}
</style>
<script>
function normaliseCellValue(value) {{
  const number = Number(value.replace(/,/g, ""));
  return Number.isNaN(number) ? value.toLowerCase() : number;
}}
function initialiseInteractiveTables() {{
  document.querySelectorAll("table[data-interactive='true']").forEach((table) => {{
    const input = document.createElement("input");
    input.className = "table-filter";
    input.type = "search";
    input.placeholder = "Filter table";
    input.setAttribute("aria-label", "Filter table");
    table.parentElement.insertBefore(input, table);

    input.addEventListener("input", () => {{
      const query = input.value.toLowerCase();
      table.querySelectorAll("tbody tr").forEach((row) => {{
        row.hidden = !row.innerText.toLowerCase().includes(query);
      }});
    }});

    table.querySelectorAll("th").forEach((header, column) => {{
      header.title = "Sort column";
      header.addEventListener("click", () => {{
        const tbody = table.tBodies[0];
        const ascending = header.dataset.sortDirection !== "asc";
        table.querySelectorAll("th").forEach((th) => delete th.dataset.sortDirection);
        header.dataset.sortDirection = ascending ? "asc" : "desc";
        Array.from(tbody.rows)
          .sort((left, right) => {{
            const a = normaliseCellValue(left.cells[column].innerText.trim());
            const b = normaliseCellValue(right.cells[column].innerText.trim());
            if (a === b) return 0;
            return (a > b ? 1 : -1) * (ascending ? 1 : -1);
          }})
          .forEach((row) => tbody.appendChild(row));
      }});
    }});
  }});
}}
function initialiseTabs() {{
  document.querySelectorAll(".tabs").forEach((tabs) => {{
    const buttons = Array.from(tabs.querySelectorAll("[role='tab']"));
    const panels = Array.from(tabs.querySelectorAll("[role='tabpanel']"));
    buttons.forEach((button, index) => {{
      button.addEventListener("click", () => {{
        buttons.forEach((item) => item.setAttribute("aria-selected", "false"));
        panels.forEach((panel) => panel.hidden = true);
        button.setAttribute("aria-selected", "true");
        panels[index].hidden = false;
        window.dispatchEvent(new Event("resize"));
      }});
    }});
  }});
}}
document.addEventListener("DOMContentLoaded", initialiseInteractiveTables);
document.addEventListener("DOMContentLoaded", initialiseTabs);
</script>
</head>
<body>
<header>
<h1>{REPORT_TITLE}</h1>
<p>This report summarises the major pre-processing and filtering steps performed by the pipeline.</p>
<nav>{nav}</nav>
</header>
<main>
{rendered_sections}
</main>
</body>
</html>
"""


def render_section(index: int, item: dict[str, str], report_text: dict[str, str]) -> str:
    title = item["title"]
    description = report_text.get(title, "")
    return f"""<section id="section-{index}">
<h2>{title}</h2>
<div class="description">{description}</div>
<div class="table-wrap">{item["body"]}</div>
</section>"""


def strip_tags(value: str) -> str:
    return value.replace("<em>", "").replace("</em>", "")


def main(snakemake):
    build_report(
        output=snakemake.output[0],
        fastq_deduplication_paths=snakemake.input.fastq_deduplication,
        fastq_trimming_path=snakemake.input.fastq_trimming,
        fastq_flash_path=snakemake.input.fastq_flash,
        fastq_digestion_paths=snakemake.input.fastq_digestion,
        reporter_filtering_paths=snakemake.input.reporters_filtering,
        reporter_deduplication_paths=snakemake.input.reporters_deduplication,
        reporter_cis_trans_paths=snakemake.input.cis_and_trans_stats,
    )


if "snakemake" in globals():
    main(globals()["snakemake"])
