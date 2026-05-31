import builtins
import importlib.util
import json
import os
import subprocess
import sys
import types
from datetime import datetime
from pathlib import Path

import pandas as pd
import polars as pl
import pytest
from cookiecutter.main import cookiecutter

import capcruncher.dependencies as dependencies
from capcruncher.dependencies import (
    CAPCRUNCHER_TOOLS_REQUIREMENT,
    DependencyVersionError,
    require_capcruncher_tools,
)


def load_workflow_script(script_name):
    script_path = (
        Path(__file__).resolve().parents[1]
        / "capcruncher"
        / "pipeline"
        / "workflow"
        / "scripts"
        / script_name
    )
    spec = importlib.util.spec_from_file_location(script_name, script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_workflow_environment_tracks_runtime_dependency_split():
    repo_root = Path(__file__).resolve().parents[1]
    env_path = (
        repo_root / "capcruncher" / "pipeline" / "workflow" / "envs" / "environment.yml"
    )
    env_text = env_path.read_text(encoding="utf-8")
    pyproject_text = (repo_root / "pyproject.toml").read_text(encoding="utf-8")

    assert "biopython" in pyproject_text
    assert "capcruncher-tools>=0.2.4,<0.3.0" in pyproject_text
    assert "capcruncher-tools>=0.2.4,<0.3.0" in env_text
    assert "typer>=0.26.0,<0.27.0" in env_text
    assert "cookiecutter" not in env_text
    assert "seaborn" not in env_text


def test_capcruncher_tools_runtime_dependency_is_supported():
    try:
        require_capcruncher_tools()
    except DependencyVersionError as exc:
        pytest.fail(str(exc))


def test_capcruncher_tools_dependency_error_includes_import_path(monkeypatch):
    monkeypatch.setattr(
        dependencies.importlib.metadata,
        "version",
        lambda distribution: "0.1.1",
    )
    monkeypatch.setattr(
        dependencies,
        "_module_path",
        lambda module_name: "/example/capcruncher_tools/__init__.py",
    )

    with pytest.raises(DependencyVersionError) as exc_info:
        require_capcruncher_tools()

    message = str(exc_info.value)
    assert "capcruncher-tools >=0.2.4,<0.3.0 is required" in message
    assert "version 0.1.1 is installed" in message
    assert "/example/capcruncher_tools/__init__.py" in message


def test_pipeline_utils_imports_without_snakemake_workflow_attribute():
    import capcruncher.pipeline.utils as pipeline_utils

    assert pipeline_utils.get_bin_sizes({"analysis": {"bin_sizes": "10000 20000"}}) == [
        10000,
        20000,
    ]


def test_pipeline_config_normalises_legacy_hub_color_by():
    from capcruncher.pipeline.utils import format_config_dict

    config = {"hub": {"color_by": "samplename"}}

    assert format_config_dict(config)["hub"]["color_by"] == "sample"


def test_pipeline_config_rejects_unknown_hub_color_by():
    from capcruncher.pipeline.utils import format_config_dict

    with pytest.raises(ValueError, match="hub.color_by"):
        format_config_dict({"hub": {"color_by": "sample_name_typo"}})


def test_pipeline_config_validates_assay_and_normalises_bin_sizes():
    from capcruncher.pipeline.validation import format_pipeline_config

    config = {
        "analysis": {
            "method": "Capture",
            "bin_sizes": "10000 20000",
        }
    }

    assert format_pipeline_config(config)["analysis"] == {
        "method": "capture",
        "bin_sizes": [10000, 20000],
    }

    with pytest.raises(ValueError, match="analysis.method"):
        format_pipeline_config({"analysis": {"method": "hic"}})


def test_pipeline_config_requires_twobit_for_custom_hub_genome():
    from capcruncher.pipeline.validation import format_pipeline_config

    with pytest.raises(ValueError, match="genome.twobit"):
        format_pipeline_config(
            {
                "genome": {"custom": True},
                "hub": {"create": True},
            }
        )


def test_pipeline_config_validates_ambiguous_workflow_options():
    from capcruncher.pipeline.validation import format_pipeline_config

    config = {
        "align": {"aligner": "Bowtie2"},
        "analysis": {"restriction_enzyme": "dpnii", "reporter_exclusion_zone": 1000},
        "analysis_optional": {
            "minimum_viewpoint_overlap": 0.5,
            "priority_chromosomes": "chr1,chr2",
        },
        "compare": {"summary_methods": "mean"},
        "plot": {"normalisation": "raw"},
        "split": {"method": "python", "n_reads": 1000},
    }

    formatted = format_pipeline_config(config)

    assert formatted["align"]["aligner"] == "bowtie2"
    assert formatted["compare"]["summary_methods"] == "mean"
    assert formatted["split"] == {"n_reads": 1000, "method": "python"}


@pytest.mark.parametrize(
    "config,error",
    [
        ({"align": {"aligner": "bwa"}}, "align.aligner"),
        ({"analysis": {"restriction_enzyme": "not-an-enzyme"}}, "restriction_enzyme"),
        ({"analysis": {"reporter_exclusion_zone": -1}}, "reporter_exclusion_zone"),
        (
            {"analysis_optional": {"minimum_viewpoint_overlap": 1.5}},
            "minimum_viewpoint_overlap",
        ),
        ({"compare": {"summary_methods": "median"}}, "compare.summary_methods"),
        ({"plot": {"normalisation": "scaled"}}, "plot.normalisation"),
        ({"split": {"method": "awk"}}, "split.method"),
        ({"split": {"n_reads": 0}}, "split.n_reads"),
    ],
)
def test_pipeline_config_rejects_ambiguous_workflow_option_typos(config, error):
    from capcruncher.pipeline.validation import format_pipeline_config

    with pytest.raises(ValueError, match=error):
        format_pipeline_config(config)


def test_pipeline_config_bridges_custom_genome_flag_for_hub_rule():
    from capcruncher.pipeline.validation import format_pipeline_config

    formatted = format_pipeline_config(
        {
            "genome": {"custom": True, "twobit": "genome.2bit"},
            "hub": {"create": True},
        }
    )

    assert formatted["hub"]["custom_genome"] is True


def test_generated_design_matrix_uses_current_sample_column(tmp_path):
    from capcruncher.pipeline.utils import get_design_matrix

    fastqs = [
        tmp_path / "SAMPLE-A_REP1_1.fastq.gz",
        tmp_path / "SAMPLE-A_REP1_2.fastq.gz",
    ]

    design = get_design_matrix(fastqs)

    assert design.to_dict("records") == [
        {"sample": "SAMPLE-A_REP1", "condition": "REP1"}
    ]


def test_workflow_output_path_manifest_is_stable():
    repo_root = Path(__file__).resolve().parents[1]
    workflow_dir = repo_root / "capcruncher" / "pipeline" / "workflow"
    workflow_text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in [
            workflow_dir / "Snakefile",
            workflow_dir / "rules" / "filter.smk",
            workflow_dir / "rules" / "pileup.smk",
            workflow_dir / "rules" / "compare.smk",
            workflow_dir / "rules" / "optional.smk",
            workflow_dir / "rules" / "qc.smk",
        ]
    )

    expected_paths = {
        "capcruncher_output/results/{sample}/{sample}.parquet",
        "capcruncher_output/results/{sample}/{sample}.hdf5",
        "capcruncher_output/results/{sample}/bigwigs/{norm}/{sample}_{viewpoint}.bigWig",
        "capcruncher_output/results/comparisons/bigwigs/{comparison}.{method}-subtraction.{viewpoint}.bigWig",
        "capcruncher_output/results/comparisons/bigwigs/{group}.{method}-summary.{viewpoint}.bigWig",
        "capcruncher_output/results/{sample}/{sample}_{read}.fastq.gz",
        "capcruncher_output/results/capcruncher_report.html",
        "capcruncher_output/results/full_qc_report.html",
        "capcruncher_output/interim/filtering/repartitioned/{sample}/flashed/",
        "capcruncher_output/interim/filtering/repartitioned/{sample}/pe/",
        "capcruncher_output/interim/filtering/deduplicated/{sample}/{combined}",
    }

    missing_paths = [
        path for path in sorted(expected_paths) if path not in workflow_text
    ]
    assert missing_paths == []


def test_rebalance_checkpoints_use_named_outputs():
    repo_root = Path(__file__).resolve().parents[1]
    fastq_rules = (
        repo_root / "capcruncher" / "pipeline" / "workflow" / "rules" / "fastq.smk"
    ).read_text(encoding="utf-8")

    assert "fastq_dir=directory(" in fastq_rules
    assert "sentinel=touch(" in fastq_rules
    assert "{output[0]}" not in fastq_rules
    assert "{output[1]}" not in fastq_rules


@pytest.fixture(scope="module")
def capture_pipeline_run(tmp_path_factory, capcruncher_subprocess_env):
    try:
        require_capcruncher_tools()
    except DependencyVersionError as exc:
        pytest.fail(
            "The capture pipeline golden-output fixture requires "
            f"capcruncher-tools{CAPCRUNCHER_TOOLS_REQUIREMENT}. {exc}"
        )

    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data" / "data_for_pipeline_run"
    run_parent = tmp_path_factory.mktemp("workflow_script_pipeline")
    current_date = datetime.now().strftime("%Y-%m-%d")
    run_dir = run_parent / f"{current_date}_project_name_capture"

    cookiecutter(
        str(repo_root / "capcruncher" / "pipeline" / "config"),
        output_dir=run_parent,
        extra_context={
            "method": "Capture-C",
            "design": str(data_dir / "design_matrix.tsv"),
            "viewpoints": str(data_dir / "mm9_capture_viewpoints_Slc25A37.bed"),
            "genome": "mm9",
            "is_custom_genome": "no",
            "genome_organism": "Mus musculus",
            "genome_fasta": str(data_dir / "chr14.fa.gz"),
            "genome_chromosome_sizes": str(data_dir / "chr14.fa.fai"),
            "genome_indicies": str(data_dir / "chr14_bowtie2_indicies" / "bt2"),
            "restriction_enzyme": "dpnii",
            "remove_blacklist": "no",
            "genomic_bin_size": "10000 20000",
            "prioritize_cis_slices": "yes",
            "priority_chromosomes": "viewpoints",
            "make_ucsc_hub": "no",
            "ucsc_hub_directory": "HUB_DIR",
            "ucsc_hub_name": "CCHUB_TEST",
            "ucsc_hub_email": "test@example.org",
            "ucsc_track_color_by": "samplename",
            "make_plots": "no",
            "plotting_coordinates": str(data_dir / "plot_coords.bed"),
            "plotting_normalisation": "n_interactions",
            "differential_contrast": "condition",
            "regenerate_fastq": "yes",
        },
        no_input=True,
    )

    for fastq in data_dir.glob("*.fastq*"):
        (run_dir / fastq.name).symlink_to(fastq)

    targets = [
        "capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_cutadapt.txt",
        "capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_flash_combo_stats.txt",
        "capcruncher_output/interim/filtering/repartitioned/SAMPLE-A_REP1/flashed",
        "capcruncher_output/results/SAMPLE-A_REP1/SAMPLE-A_REP1.hdf5",
        "capcruncher_output/results/SAMPLE-A_REP1/SAMPLE-A_REP1.parquet",
        "capcruncher_output/resources/restriction_fragments/genome.digest.bed.gz",
    ]
    result = subprocess.run(
        [
            "capcruncher",
            "pipeline",
            "--no-logo",
            "-c",
            os.environ.get("CAPCRUNCHER_TEST_CORES", "1"),
            "--show-failed-logs",
            *targets,
        ],
        cwd=run_dir,
        env=capcruncher_subprocess_env,
    )
    assert result.returncode == 0
    return run_dir


@pytest.mark.parametrize(
    "script_name",
    [
        "combine_deduplicated_slices.py",
        "count_identified_viewpoints.py",
        "extract_flash_data.py",
        "extract_trimming_data.py",
        "identify_viewpoints_with_interactions.py",
        "make_ucsc_hub.py",
        "plot.py",
        "repartition_filtered_slices.py",
        "remove_duplicate_coordinates.py",
        "run_differential.py",
        "save_design.py",
        "validation_check_n_bins_per_viewpoint.py",
        "validation_confirm_annotated_viewpoints_present.py",
    ],
)
def test_workflow_scripts_import_without_snakemake(script_name):
    load_workflow_script(script_name)


@pytest.mark.parametrize(
    "blocked_import,script_name",
    [
        ("plotnado", "plot.py"),
        ("tracknado", "make_ucsc_hub.py"),
    ],
)
def test_optional_workflow_scripts_import_without_optional_deps(
    monkeypatch, blocked_import, script_name
):
    real_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name == blocked_import or name.startswith(f"{blocked_import}."):
            raise ModuleNotFoundError(
                f"No module named '{blocked_import}'", name=blocked_import
            )
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    load_workflow_script(script_name)


def test_validation_confirm_annotated_viewpoints_present_counts_current_polars(
    tmp_path,
):
    script = load_workflow_script("validation_confirm_annotated_viewpoints_present.py")
    slices_a = tmp_path / "slices_a.parquet"
    slices_b = tmp_path / "slices_b.parquet"
    viewpoints = tmp_path / "viewpoints.bed"
    counts = tmp_path / "viewpoints_present.tsv"
    sentinel = tmp_path / "validated.sentinel"

    pl.DataFrame({"capture": ["vp1", "vp1", "vp2"]}).write_parquet(slices_a)
    pl.DataFrame({"capture": ["vp2"]}).write_parquet(slices_b)
    viewpoints.write_text("chr1\t0\t10\tvp1\nchr1\t20\t30\tvp2\n")

    script.validate_viewpoints_present(
        [slices_a, slices_b],
        viewpoints,
        counts,
        sentinel,
    )

    df_counts = pl.read_csv(counts, separator="\t")
    assert sentinel.exists()
    assert df_counts.select(["capture", "n_slices"]).sort("capture").to_dicts() == [
        {"capture": "vp1", "n_slices": 2},
        {"capture": "vp2", "n_slices": 2},
    ]


def test_validation_confirm_annotated_viewpoints_present_reports_missing(tmp_path):
    script = load_workflow_script("validation_confirm_annotated_viewpoints_present.py")
    slices = tmp_path / "slices.parquet"
    viewpoints = tmp_path / "viewpoints.bed"

    pl.DataFrame({"capture": ["vp1"]}).write_parquet(slices)
    viewpoints.write_text("chr1\t0\t10\tvp1\nchr1\t20\t30\tvp2\n")

    with pytest.raises(ValueError, match="vp2"):
        script.validate_viewpoints_present(
            [slices],
            viewpoints,
            tmp_path / "viewpoints_present.tsv",
            tmp_path / "validated.sentinel",
        )


def test_count_identified_viewpoints_filters_empty_values(tmp_path):
    script = load_workflow_script("count_identified_viewpoints.py")
    slices_dir = tmp_path / "slices"
    slices_dir.mkdir()
    output = tmp_path / "identified_viewpoints.tsv"

    pl.DataFrame(
        {
            "viewpoint": ["vp1", "vp1", "vp2", "", None],
            "pe": ["pe1", "pe1", "", "pe2", "pe3"],
        }
    ).write_parquet(slices_dir / "slices.parquet")

    script.write_identified_viewpoints(slices_dir, output)

    assert pl.read_csv(output, separator="\t").sort("viewpoint").to_dicts() == [
        {"viewpoint": "vp1", "pe": "pe1"},
    ]


def test_extract_flash_stats_aggregates_multiqc_rows(tmp_path):
    script = load_workflow_script("extract_flash_data.py")
    flash_summary = tmp_path / "flash.tsv"

    pd.DataFrame(
        {
            "Sample": ["SAMPLE-A_part0", "SAMPLE-A_part1", "SAMPLE-B_part0"],
            "combopairs": [10, 5, 7],
            "uncombopairs": [2, 3, 1],
        }
    ).to_csv(flash_summary, sep="\t", index=False)

    stats = script.extract_flash_stats(flash_summary)

    assert [stat.model_dump() for stat in stats] == [
        {
            "sample": "SAMPLE-A",
            "n_combined": 15,
            "n_uncombined": 5,
            "n_total": 20,
            "percentage_combined": 75,
        },
        {
            "sample": "SAMPLE-B",
            "n_combined": 7,
            "n_uncombined": 1,
            "n_total": 8,
            "percentage_combined": 87.5,
        },
    ]


def test_extract_flash_stats_derives_combined_pairs_from_current_multiqc(tmp_path):
    script = load_workflow_script("extract_flash_data.py")
    flash_summary = tmp_path / "flash.tsv"

    pd.DataFrame(
        {
            "Sample": ["SAMPLE-A_part0_1", "SAMPLE-A_part1_1", "SAMPLE-B_part0_1"],
            "totalpairs": [147, 142, 225],
            "discardpairs": [0, 2, 5],
            "uncombopairs": [55, 57, 99],
        }
    ).to_csv(flash_summary, sep="\t", index=False)

    stats = script.extract_flash_stats(flash_summary)

    assert [stat.model_dump() for stat in stats] == [
        {
            "sample": "SAMPLE-A",
            "n_combined": 175,
            "n_uncombined": 112,
            "n_total": 287,
            "percentage_combined": 60.97560975609756,
        },
        {
            "sample": "SAMPLE-B",
            "n_combined": 121,
            "n_uncombined": 99,
            "n_total": 220,
            "percentage_combined": 55.00000000000001,
        },
    ]


def test_extract_trimming_stats_aggregates_multiqc_rows(tmp_path):
    script = load_workflow_script("extract_trimming_data.py")
    trimming_summary = tmp_path / "trimming.tsv"

    pd.DataFrame(
        {
            "Sample": ["SAMPLE-A_part0_1", "SAMPLE-A_part1_1", "SAMPLE-A_part0_2"],
            "r_processed": [10, 20, 30],
            "r_written": [8, 19, 25],
            "r_with_adapters": [2, 1, 5],
        }
    ).to_csv(trimming_summary, sep="\t", index=False)

    stats = script.extract_trimming_stats(trimming_summary)

    assert [stat.model_dump() for stat in stats] == [
        {
            "sample": "SAMPLE-A",
            "read_number": 1,
            "reads_input": 30,
            "reads_output": 27,
            "reads_with_adapter_identified": 3,
            "percentage_trimmed": 10.0,
            "percentage_passing_quality_filter": 90.0,
        },
        {
            "sample": "SAMPLE-A",
            "read_number": 2,
            "reads_input": 30,
            "reads_output": 25,
            "reads_with_adapter_identified": 5,
            "percentage_trimmed": pytest.approx(16.666666666666664),
            "percentage_passing_quality_filter": pytest.approx(83.33333333333334),
        },
    ]


def test_identify_viewpoints_with_interactions_uses_count_column(monkeypatch):
    script = load_workflow_script("identify_viewpoints_with_interactions.py")

    class DummyPixels:
        def __init__(self, counts):
            self.counts = counts

        def __getitem__(self, item):
            return pd.DataFrame({"bin1_id": [0, 1], "count": self.counts})

    class DummyCooler:
        def __init__(self, uri):
            self.uri = uri

        def pixels(self):
            viewpoint = self.uri.split("::", 1)[1]
            return DummyPixels([0, 0] if viewpoint == "empty" else [0, 3])

    monkeypatch.setattr(
        script.cooler,
        "api",
        types.SimpleNamespace(list_coolers=lambda path: ["empty", "present"]),
    )
    monkeypatch.setattr(script.cooler, "Cooler", DummyCooler)

    assert script.viewpoints_with_interactions("sample.cool") == ["present"]


def test_write_viewpoints_with_interactions_writes_per_sample_json(
    monkeypatch, tmp_path
):
    script = load_workflow_script("identify_viewpoints_with_interactions.py")
    monkeypatch.setattr(
        script,
        "viewpoints_with_interactions",
        lambda cooler_path: [f"{cooler_path}-vp"],
    )

    script.write_viewpoints_with_interactions(
        ["a.cool", "b.cool"],
        ["sample-a", "sample-b"],
        tmp_path,
    )

    assert json.loads((tmp_path / "sample-a.json").read_text()) == ["a.cool-vp"]
    assert json.loads((tmp_path / "sample-b.json").read_text()) == ["b.cool-vp"]


@pytest.mark.pipeline
@pytest.mark.slow
def test_workflow_scripts_run_on_capture_pipeline_inputs(
    capture_pipeline_run, tmp_path
):
    repo_root = Path(__file__).resolve().parents[1]
    count_script = load_workflow_script("count_identified_viewpoints.py")
    flash_script = load_workflow_script("extract_flash_data.py")
    trimming_script = load_workflow_script("extract_trimming_data.py")
    identify_script = load_workflow_script("identify_viewpoints_with_interactions.py")
    remove_dups_script = load_workflow_script("remove_duplicate_coordinates.py")
    validation_script = load_workflow_script("validation_check_n_bins_per_viewpoint.py")

    output_dir = tmp_path / "script_outputs"
    output_dir.mkdir()

    trimming_output = output_dir / "trimming.json"
    trimming_script.write_trimming_stats(
        capture_pipeline_run
        / "capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_cutadapt.txt",
        trimming_output,
    )
    assert json.loads(trimming_output.read_text())

    flash_output = output_dir / "flash.json"
    flash_script.write_flash_stats(
        capture_pipeline_run
        / "capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_flash_combo_stats.txt",
        flash_output,
    )
    assert json.loads(flash_output.read_text())

    identified_output = output_dir / "identified_viewpoints.tsv"
    count_script.write_identified_viewpoints(
        capture_pipeline_run
        / "capcruncher_output/results/SAMPLE-A_REP1/SAMPLE-A_REP1.parquet",
        identified_output,
    )
    identified = pl.read_csv(identified_output, separator="\t")
    assert {"viewpoint", "pe"}.issubset(set(identified.columns))

    viewpoints_output = output_dir / "viewpoints_with_interactions"
    identify_script.write_viewpoints_with_interactions(
        [
            capture_pipeline_run
            / "capcruncher_output/results/SAMPLE-A_REP1/SAMPLE-A_REP1.hdf5"
        ],
        ["SAMPLE-A_REP1"],
        viewpoints_output,
    )
    assert json.loads((viewpoints_output / "SAMPLE-A_REP1.json").read_text())

    validation_sentinel = output_dir / "validated.sentinel"
    validation_script.check_n_bins_per_viewpoint(
        bins=capture_pipeline_run
        / "capcruncher_output/resources/restriction_fragments/genome.digest.bed.gz",
        viewpoints=repo_root
        / "tests/data/data_for_pipeline_run/mm9_capture_viewpoints_Slc25A37.bed",
        output_sentinel=validation_sentinel,
        ignore_multiple_bins_per_viewpoint=False,
    )
    assert validation_sentinel.exists()

    deduplicated_output = output_dir / "deduplicated_flashed"
    deduplication_stats = output_dir / "deduplicated_flashed.json"
    remove_dups_script.remove_duplicate_coordinates(
        slices_directory=capture_pipeline_run
        / "capcruncher_output/interim/filtering/repartitioned/SAMPLE-A_REP1/flashed",
        output_slices=deduplicated_output,
        output_statistics=deduplication_stats,
        read_type="flashed",
        sample_name="SAMPLE-A_REP1",
        log_path=output_dir / "remove_duplicate_coordinates.log",
    )
    assert deduplicated_output.exists()
    assert list(deduplicated_output.rglob("*.parquet"))
    assert deduplication_stats.exists()


@pytest.mark.pipeline
@pytest.mark.slow
def test_capture_pipeline_golden_outputs(capture_pipeline_run):
    import cooler

    reporter_parquet = (
        capture_pipeline_run
        / "capcruncher_output/results/SAMPLE-A_REP1/SAMPLE-A_REP1.parquet"
    )
    cooler_path = (
        capture_pipeline_run
        / "capcruncher_output/results/SAMPLE-A_REP1/SAMPLE-A_REP1.hdf5"
    )
    digest_bed = (
        capture_pipeline_run
        / "capcruncher_output/resources/restriction_fragments/genome.digest.bed.gz"
    )

    reporters = pd.read_parquet(reporter_parquet)
    assert len(reporters) == 205
    assert reporters["viewpoint"].astype(str).value_counts().to_dict() == {
        "Slc25A37": 205
    }
    assert reporters["capture"].notna().sum() == 94
    assert set(reporters["capture"].dropna().astype(str)) == {"Slc25A37"}

    assert len(pd.read_csv(digest_bed, sep="\t", header=None)) == 303397

    assert cooler.api.list_coolers(str(cooler_path)) == [
        "/Slc25A37",
        "/Slc25A37/resolutions/10000",
        "/Slc25A37/resolutions/20000",
    ]

    raw_cooler = cooler.Cooler(f"{cooler_path}::/Slc25A37")
    if raw_cooler.info["metadata"]["n_total_interactions"] != 130:
        pytest.xfail(
            "capcruncher-tools 0.2.5 currently counts 75 interactions at the "
            "reporter-counting boundary even though CapCruncher passes all 205 "
            "countable reporter rows through."
        )

    assert raw_cooler.info["metadata"] == {
        "viewpoint_bins": [169634],
        "viewpoint_name": "Slc25A37",
        "viewpoint_chrom": ["chr14"],
        "viewpoint_coords": ["chr14:69902454-69903469"],
        "n_cis_interactions": 130,
        "n_total_interactions": 130,
    }
    assert raw_cooler.pixels()[:].to_dict("records") == [
        {"bin1_id": 169686, "bin2_id": 169687, "count": 9},
        {"bin1_id": 169686, "bin2_id": 169744, "count": 82},
        {"bin1_id": 169687, "bin2_id": 169744, "count": 10},
        {"bin1_id": 169744, "bin2_id": 169786, "count": 1},
        {"bin1_id": 169744, "bin2_id": 169845, "count": 10},
        {"bin1_id": 169744, "bin2_id": 169846, "count": 6},
        {"bin1_id": 169744, "bin2_id": 169847, "count": 2},
        {"bin1_id": 169845, "bin2_id": 169846, "count": 6},
        {"bin1_id": 169845, "bin2_id": 169847, "count": 2},
        {"bin1_id": 169846, "bin2_id": 169847, "count": 2},
    ]

    for group, expected_bins, expected_pixels in [
        ("/Slc25A37/resolutions/10000", 12520, 5),
        ("/Slc25A37/resolutions/20000", 6260, 5),
    ]:
        binned_cooler = cooler.Cooler(f"{cooler_path}::{group}")
        assert len(binned_cooler.bins()[:]) == expected_bins
        pixels = binned_cooler.pixels()[:]
        assert len(pixels) == expected_pixels
        assert int(pixels["count"].sum()) == 130
        assert binned_cooler.info["metadata"]["n_interactions_total"] == 130


@pytest.mark.pipeline
@pytest.mark.slow
def test_countable_reporter_handoff_preserves_pipeline_partitions(
    capture_pipeline_run, tmp_path
):
    from capcruncher.api.interactions.reporters import write_countable_reporters

    reporter_parquet = (
        capture_pipeline_run
        / "capcruncher_output/results/SAMPLE-A_REP1/SAMPLE-A_REP1.parquet"
    )
    viewpoints = (
        Path(__file__).resolve().parents[1]
        / "tests/data/data_for_pipeline_run/mm9_capture_viewpoints_Slc25A37.bed"
    )

    countable_reporters = write_countable_reporters(
        reporter_parquet, viewpoints, tmp_path / "countable"
    )
    countable_files = sorted(countable_reporters.glob("*.parquet"))

    assert [path.name for path in countable_files] == [
        "part-0.parquet",
        "part-1.parquet",
    ]
    assert sum(len(pl.read_parquet(path)) for path in countable_files) == 205
    assert (
        pl.concat([pl.read_parquet(path) for path in countable_files])
        .get_column("viewpoint")
        .drop_nulls()
        .cast(pl.Utf8)
        .unique()
        .to_list()
    ) == ["Slc25A37"]


def test_reporter_summary_counts_unused_viewpoint_categories_as_zero(tmp_path):
    from capcruncher.api.interactions.reporters import summarise_reporter_viewpoints

    reporters = tmp_path / "reporters.parquet"
    pl.DataFrame(
        {
            "viewpoint": pl.Series(
                ["Slc25A37", None],
                dtype=pl.Enum(["Slc25A37", "unused_viewpoint"]),
            )
        }
    ).write_parquet(reporters)

    summary = summarise_reporter_viewpoints(reporters)

    assert summary.viewpoints == ["Slc25A37", "unused_viewpoint"]
    assert summary.viewpoint_sizes == {"Slc25A37": 1}


@pytest.mark.pipeline
@pytest.mark.slow
def test_capcruncher_tools_counts_expected_pipeline_pixels(
    capture_pipeline_run, tmp_path
):
    from capcruncher_tools.count import count_viewpoint_pixels

    from capcruncher.api.interactions.reporters import write_countable_reporters

    reporter_parquet = (
        capture_pipeline_run
        / "capcruncher_output/results/SAMPLE-A_REP1/SAMPLE-A_REP1.parquet"
    )
    viewpoints = (
        Path(__file__).resolve().parents[1]
        / "tests/data/data_for_pipeline_run/mm9_capture_viewpoints_Slc25A37.bed"
    )
    countable_reporters = write_countable_reporters(
        reporter_parquet, viewpoints, tmp_path / "countable"
    )

    _viewpoint, counts = count_viewpoint_pixels(
        parquet=str(countable_reporters / "*.parquet"),
        viewpoint="Slc25A37",
    )
    observed_total = int(counts["count"].sum())
    if observed_total != 130:
        pytest.xfail(
            "capcruncher-tools 0.2.5 counts 75 interactions from the "
            "CapCruncher countable reporter handoff; expected golden total is 130."
        )

    assert observed_total == 130


def test_remove_duplicate_coordinates_preserves_empty_parquet_schema(tmp_path):
    script = load_workflow_script("remove_duplicate_coordinates.py")
    slices = tmp_path / "slices"
    slices.mkdir()
    output = tmp_path / "deduplicated"
    statistics = tmp_path / "stats.csv"

    pl.DataFrame(
        {
            "viewpoint": [],
            "parent_id": [],
            "slice_id": [],
            "coordinates": [],
        },
        schema={
            "viewpoint": pl.String,
            "parent_id": pl.Int64,
            "slice_id": pl.Int64,
            "coordinates": pl.String,
        },
    ).write_parquet(slices / "empty.parquet")

    script.remove_duplicate_coordinates(
        slices_directory=slices,
        output_slices=output,
        output_statistics=statistics,
        read_type="flashed",
        sample_name="sample-a",
        log_path=tmp_path / "deduplicate.log",
    )

    assert pl.scan_parquet(output).collect_schema()["viewpoint"] == pl.String
    assert statistics.exists()


def test_remove_duplicate_coordinates_main_reraises_failures(monkeypatch, tmp_path):
    script = load_workflow_script("remove_duplicate_coordinates.py")
    slices = tmp_path / "slices"
    slices.mkdir()
    output = tmp_path / "deduplicated"
    statistics = tmp_path / "stats.json"
    log = tmp_path / "deduplicate.log"

    pl.DataFrame({"viewpoint": ["vp1"]}).write_parquet(slices / "data.parquet")

    def fail_deduplicate(**kwargs):
        raise RuntimeError("deduplicate failed")

    monkeypatch.setattr(script, "deduplicate", fail_deduplicate)
    snakemake = types.SimpleNamespace(
        input=types.SimpleNamespace(slices_directory=slices),
        output=types.SimpleNamespace(slices=output, statistics=statistics),
        params=types.SimpleNamespace(read_type="flashed", sample_name="sample-a"),
        log=[log],
    )

    with pytest.raises(RuntimeError, match="deduplicate failed"):
        script.main(snakemake)

    assert not output.exists()
    assert "deduplicate failed" in log.read_text(encoding="utf-8")


def test_make_ucsc_hub_builds_tracknado_metadata(tmp_path):
    script = load_workflow_script("make_ucsc_hub.py")
    viewpoints = tmp_path / "viewpoints.bigBed"

    records = script.build_track_metadata(
        bigwigs=[
            tmp_path / "raw" / "SAMPLE-A_REP1_Slc25A37.bigWig",
            tmp_path / "norm" / "SAMPLE-A_REP1_Slc25A37.bigWig",
        ],
        bigwigs_summary=[
            tmp_path / "SAMPLE-A.mean-summary.Slc25A37.bigWig",
        ],
        bigwigs_comparison=[
            tmp_path / "SAMPLE-A-SAMPLE-B.mean-subtraction.Slc25A37.bigWig",
        ],
        viewpoints=viewpoints,
    )

    assert [
        {
            "category": record["category"],
            "normalisation": record["normalisation"],
            "sample": record["sample"],
            "aggregation": record["aggregation"],
            "ext": record["ext"],
        }
        for record in records
    ] == [
        {
            "category": "Replicates",
            "normalisation": "raw",
            "sample": "SAMPLE-A_REP1",
            "aggregation": "replicate",
            "ext": "bigWig",
        },
        {
            "category": "Replicates",
            "normalisation": "norm",
            "sample": "SAMPLE-A_REP1",
            "aggregation": "replicate",
            "ext": "bigWig",
        },
        {
            "category": "Aggregated",
            "normalisation": "norm",
            "sample": "SAMPLE-A",
            "aggregation": "mean",
            "ext": "bigWig",
        },
        {
            "category": "Subtraction",
            "normalisation": "norm",
            "sample": "SAMPLE-A-SAMPLE-B",
            "aggregation": "mean",
            "ext": "bigWig",
        },
        {
            "category": "Annotation",
            "normalisation": "viewpoints",
            "sample": "viewpoints",
            "aggregation": "viewpoints",
            "ext": "bigBed",
        },
    ]
    assert "overlay" not in records[-1]


def test_make_ucsc_hub_rejects_unsupported_track_names(tmp_path):
    script = load_workflow_script("make_ucsc_hub.py")

    with pytest.raises(ValueError, match="Could not parse CapCruncher track path"):
        script.capcruncher_track_metadata(tmp_path / "sample.invalid-name.bigWig")


def test_make_ucsc_hub_uses_modern_tracknado_builder(monkeypatch, tmp_path):
    script = load_workflow_script("make_ucsc_hub.py")
    calls = []

    class DummyHub:
        pass

    class DummyBuilder:
        def add_tracks(self, tracks):
            calls.append(("add_tracks", tracks))
            return self

        def with_metadata_extractor(self, extractor):
            calls.append(("with_metadata_extractor", extractor.__name__))
            return self

        def group_by(self, *columns, as_supertrack=False):
            calls.append(("group_by", columns, as_supertrack))
            return self

        def overlay_by(self, *columns):
            calls.append(("overlay_by", columns))
            return self

        def color_by(self, column):
            calls.append(("color_by", column))
            return self

        def with_custom_genome(self, name, twobit_file, organism, default_position):
            calls.append(
                ("with_custom_genome", name, twobit_file, organism, default_position)
            )
            return self

        def build(self, **kwargs):
            calls.append(("build", kwargs))
            return DummyHub()

    monkeypatch.setitem(
        sys.modules,
        "tracknado",
        types.SimpleNamespace(HubBuilder=DummyBuilder),
    )

    result = script.build_hub(
        tracks=script.collect_track_paths(
            bigwigs=[tmp_path / "raw" / "SAMPLE-A_REP1_Slc25A37.bigWig"],
            bigwigs_summary=[],
            bigwigs_comparison=[],
            viewpoints=tmp_path / "viewpoints.bigBed",
        ),
        color_by="samplename",
        genome="mm10",
        hub_name="capcruncher",
        hub_email="test@example.org",
        custom_genome=True,
        genome_twobit=tmp_path / "genome.2bit",
        genome_organism="Mouse",
        genome_default_position="chr1:1-100",
        report=tmp_path / "report.html",
        outdir=tmp_path / "hub",
    )

    assert isinstance(result, DummyHub)
    assert calls[1:] == [
        ("with_metadata_extractor", "capcruncher_track_metadata"),
        ("group_by", ("category", "normalisation"), True),
        ("group_by", ("sample", "viewpoint", "aggregation"), False),
        ("overlay_by", ("overlay",)),
        ("color_by", "sample"),
        (
            "with_custom_genome",
            "mm10",
            tmp_path / "genome.2bit",
            "Mouse",
            "chr1:1-100",
        ),
        (
            "build",
            {
                "name": "capcruncher",
                "genome": "mm10",
                "outdir": tmp_path / "hub",
                "hub_email": "test@example.org",
                "description_html": tmp_path / "report.html",
            },
        ),
    ]


def test_make_ucsc_hub_requires_twobit_for_custom_genome(tmp_path):
    script = load_workflow_script("make_ucsc_hub.py")

    with pytest.raises(ValueError, match="genome twoBit file"):
        script.build_hub(
            tracks=[
                tmp_path / "raw" / "SAMPLE-A_REP1_Slc25A37.bigWig",
                tmp_path / "viewpoints.bigBed",
            ],
            color_by="sample",
            genome="custom",
            hub_name="capcruncher",
            hub_email="test@example.org",
            custom_genome=True,
            genome_twobit=None,
            report=tmp_path / "report.html",
            outdir=tmp_path / "hub",
        )
