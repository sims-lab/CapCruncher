import json

from capcruncher.pipeline.workflow.report.make_report import build_report


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_report_builds_interactive_html_without_quarto(tmp_path):
    stats = tmp_path / "stats"
    output = tmp_path / "capcruncher_report.html"

    fastq_deduplication = [
        write_json(
            stats / "deduplication" / "data" / "SAMPLE.deduplication.json",
            {"sample": "SAMPLE", "total": 100, "duplicates": 10},
        )
    ]
    trimming = write_json(
        stats / "trimming" / "trimming.json",
        [
            json.dumps(
                {
                    "sample": "SAMPLE",
                    "read_number": 1,
                    "reads_input": 90,
                    "reads_output": 85,
                    "reads_with_adapter_identified": 5,
                }
            ),
            json.dumps(
                {
                    "sample": "SAMPLE",
                    "read_number": 2,
                    "reads_input": 90,
                    "reads_output": 84,
                    "reads_with_adapter_identified": 6,
                }
            ),
        ],
    )
    flash = write_json(
        stats / "flash" / "flash.json",
        [json.dumps({"sample": "SAMPLE", "n_combined": 60, "n_uncombined": 30})],
    )
    digestion = [
        write_json(
            stats / "digestion" / "data" / "SAMPLE_part0_flashed.json",
            {
                "sample": "SAMPLE",
                "read_type": "Flashed",
                "read_stats": {
                    "unfiltered": {"read1": 60, "read2": 0},
                    "filtered": {"read1": 55, "read2": 0},
                },
                "histograms": {
                    "lengths": {
                        "read1": {"name": "slice_length", "hist": {"25": 3}},
                        "read2": {"name": "slice_length", "hist": {}},
                    }
                },
            },
        ),
        write_json(
            stats / "digestion" / "data" / "SAMPLE_part0_pe.json",
            {
                "sample": "SAMPLE",
                "read_type": "Pe",
                "read_stats": {
                    "unfiltered": {"read1": 30, "read2": 30},
                    "filtered": {"read1": 28, "read2": 27},
                },
                "histograms": {
                    "lengths": {
                        "read1": {"name": "slice_length", "hist": {"30": 2}},
                        "read2": {"name": "slice_length", "hist": {"35": 2}},
                    }
                },
            },
        ),
    ]
    filtering = [
        write_json(
            stats / "filtering" / "data" / "SAMPLE_part0_flashed.json",
            {
                "stats": [
                    {
                        "sample": "SAMPLE",
                        "stage": "pre-filtering",
                        "n_fragments": 55,
                        "n_slices": 70,
                        "read_type": "flashed",
                    },
                    {
                        "sample": "SAMPLE",
                        "stage": "contains_capture_and_reporter",
                        "n_fragments": 45,
                        "n_slices": 60,
                        "read_type": "flashed",
                    },
                ]
            },
        )
    ]
    reporter_deduplication = [
        write_json(
            stats / "deduplication_final" / "data" / "SAMPLE_flashed.json",
            {
                "sample": "SAMPLE",
                "read_type": "flashed",
                "n_total_reads": 45,
                "n_unique_reads": 40,
                "n_total_slices": 60,
                "n_unique_slices": 52,
            },
        )
    ]
    cis_trans = [
        write_json(
            stats / "cis_and_trans_reporters" / "data" / "SAMPLE.json",
            {
                "stats": [
                    {
                        "sample": "SAMPLE",
                        "read_type": "flashed",
                        "viewpoint": "VP1",
                        "cis_or_trans": "cis",
                        "count": 40,
                    }
                ]
            },
        )
    ]

    build_report(
        output,
        fastq_deduplication,
        trimming,
        flash,
        digestion,
        filtering,
        reporter_deduplication,
        cis_trans,
    )

    html = output.read_text(encoding="utf-8")
    assert "CapCruncher Run Report" in html
    assert "Plotly.newPlot" in html
    assert 'role="tab"' in html
    assert "Filter table" in html
    assert "quarto" not in html.lower()
    assert ".qmd" not in html
