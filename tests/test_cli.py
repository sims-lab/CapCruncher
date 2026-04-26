from loguru import logger
import pandas as pd
import pytest
import os
import pathlib
import subprocess
import sys
from click.testing import CliRunner
import glob
from types import SimpleNamespace

from capcruncher.cli import cli
from capcruncher.cli import pipeline as cli_pipeline
from capcruncher.cli.interactions import _write_countable_reporters


@pytest.fixture(scope="module", autouse=True)
def setup_testing_dir(tmpdir_factory):
    cwd = os.getcwd()
    tmpdir = tmpdir_factory.mktemp("cli_testing")
    os.chdir(tmpdir)
    yield
    os.chdir(cwd)


@pytest.fixture(scope="session")
def testdata_dirname():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    return dirname


@pytest.fixture
def data_digestion(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "fastq_digestion")
    return data_dir


@pytest.fixture
def data_deduplication(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "fastq_deduplication")
    return data_dir


@pytest.fixture(scope="module")
def data_annotation(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "alignment_annotation")
    return data_dir


@pytest.fixture(scope="module")
def data_filter(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "alignment_filtering")
    return data_dir


@pytest.fixture(scope="module")
def data_deduplication_alignments(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "alignment_deduplication")
    return data_dir


@pytest.fixture(scope="module")
def data_reporters_count(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "reporters_count")
    return data_dir


@pytest.fixture(scope="module")
def data_reporters_store(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "reporters_store")
    return data_dir


@pytest.fixture(scope="module")
def data_pipeline(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "data_for_pipeline_run")
    return data_dir


@pytest.fixture(scope="module")
def cli_runner():
    return CliRunner()


def test_cli_runs(cli_runner):
    """Test checks that the cli is functional and the help option works"""

    import capcruncher.cli

    assert capcruncher.cli.cli is cli
    result = cli_runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "pipeline-init" in result.output


def test_cli_import_does_not_import_heavy_runtime_modules():
    code = (
        "import sys; "
        "import capcruncher.cli; "
        "blocked = [name for name in ('pandas', 'polars', 'pyranges1') if name in sys.modules]; "
        "print(','.join(blocked))"
    )

    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        check=True,
        text=True,
    )

    assert result.stdout.strip() == ""


@pytest.mark.parametrize(
    "command",
    [
        ["interactions", "differential", "--help"],
        ["interactions", "compare", "differential", "--help"],
    ],
)
def test_differential_help_does_not_import_pydeseq2(
    cli_runner, monkeypatch, command
):
    real_import = __import__

    def guarded_import(name, *args, **kwargs):
        if name == "pydeseq2" or name.startswith("pydeseq2."):
            raise ModuleNotFoundError("No module named 'pydeseq2'", name="pydeseq2")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr("builtins.__import__", guarded_import)

    result = cli_runner.invoke(cli, command)

    assert result.exit_code == 0


@pytest.mark.parametrize(
    "command",
    [
        ["plot", "-r", "chr1:1-100", "-t", "template.toml", "-o", "plot.png"],
        ["plot", "render", "-r", "chr1:1-100", "-t", "template.toml", "-o", "plot.png"],
    ],
)
def test_plot_render_commands(cli_runner, monkeypatch, command):
    import capcruncher.cli.plot as cli_plot

    calls = []

    def fake_render_plot(region, template, output):
        calls.append((region, template, output))

    monkeypatch.setattr(cli_plot, "render_plot", fake_render_plot)

    result = cli_runner.invoke(cli, command)

    assert result.exit_code == 0
    assert calls == [("chr1:1-100", "template.toml", "plot.png")]


def test_plot_help_does_not_import_plotnado(cli_runner, monkeypatch):
    real_import = __import__

    def guarded_import(name, *args, **kwargs):
        if name == "plotnado" or name.startswith("plotnado."):
            raise ModuleNotFoundError("No module named 'plotnado'", name="plotnado")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr("builtins.__import__", guarded_import)

    result = cli_runner.invoke(cli, ["plot", "render", "--help"])

    assert result.exit_code == 0


@pytest.mark.parametrize(
    "command",
    [
        [
            "genome",
            "digest",
            "genome.fa",
            "--recognition-site",
            "dpnii",
            "--output-file",
            "digest.bed",
            "--sort",
        ],
        [
            "genome",
            "digest",
            "genome.fa",
            "--recognition_site",
            "dpnii",
            "--output_file",
            "digest.bed",
            "--sort",
        ],
    ],
)
def test_genome_digest_option_aliases(cli_runner, monkeypatch, command):
    import capcruncher.api.genome as genome_api

    calls = []

    def fake_digest(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr(genome_api, "digest_genome", fake_digest)

    result = cli_runner.invoke(cli, command)

    assert result.exit_code == 0
    assert calls == [
        {
            "input_fasta": "genome.fa",
            "recognition_site": "dpnii",
            "output_file": "digest.bed",
            "logfile": "genome_digest.log",
            "remove_cutsite": True,
            "sort": True,
        }
    ]


def test_alignments_typer_option_aliases(cli_runner, monkeypatch):
    import capcruncher.api.alignments_annotate as alignments_annotate
    import capcruncher.api.alignments_filter as alignments_filter

    annotate_calls = []
    filter_calls = []

    def fake_annotate(**kwargs):
        annotate_calls.append(kwargs)

    def fake_filter(**kwargs):
        filter_calls.append(kwargs)

    monkeypatch.setattr(alignments_annotate, "annotate", fake_annotate)
    monkeypatch.setattr(alignments_filter, "filter", fake_filter)

    annotate_result = cli_runner.invoke(
        cli,
        [
            "alignments",
            "annotate",
            "slices.bam",
            "--bed_files",
            "targets.bed",
            "--actions",
            "get",
            "--names",
            "targets",
            "--overlap_fractions",
            "0.5",
            "--n_cores",
            "2",
            "--invalid_bed_action",
            "ignore",
        ],
    )
    filter_result = cli_runner.invoke(
        cli,
        [
            "alignments",
            "filter",
            "capture",
            "--bam",
            "reads.bam",
            "--annotations",
            "annotations.parquet",
            "--output_prefix",
            "filtered",
            "--no-fragments",
        ],
    )

    assert annotate_result.exit_code == 0
    assert filter_result.exit_code == 0
    assert annotate_calls == [
        {
            "slices": "slices.bam",
            "actions": ("get",),
            "bed_files": ("targets.bed",),
            "names": ("targets",),
            "overlap_fractions": (0.5,),
            "dtypes": ("str",),
            "output": "annotated.slices.parquet",
            "duplicates": "remove",
            "n_cores": 2,
            "invalid_bed_action": "ignore",
            "blacklist": "",
            "prioritize_cis_slices": False,
            "priority_chroms": "",
        }
    ]
    assert filter_calls == [
        {
            "method": "capture",
            "bam": "reads.bam",
            "annotations": "annotations.parquet",
            "custom_filtering": None,
            "output_prefix": "filtered",
            "statistics": "filtering_stats.json",
            "sample_name": None,
            "read_type": "flashed",
            "fragments": False,
        }
    ]


@pytest.mark.parametrize(
    "command",
    [
        ["pipeline", "config", "--help"],
        ["pipeline-config", "--help"],
    ],
)
def test_pipeline_config_help_does_not_import_cookiecutter(
    cli_runner, monkeypatch, command
):
    real_import = __import__

    def guarded_import(name, *args, **kwargs):
        if name == "cookiecutter" or name.startswith("cookiecutter."):
            raise ModuleNotFoundError(
                "No module named 'cookiecutter'", name="cookiecutter"
            )
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr("builtins.__import__", guarded_import)

    result = cli_runner.invoke(cli, command)

    assert result.exit_code == 0


def test_pipeline_init_installs_presets(cli_runner, tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))

    result = cli_runner.invoke(cli, ["pipeline-init"])

    assert result.exit_code == 0
    profiles_dir = tmp_path / "snakemake"
    assert (profiles_dir / "capcruncher-local" / "profile.v9+.yaml").exists()
    assert (profiles_dir / "capcruncher-local-conda" / "profile.v9+.yaml").exists()
    assert (profiles_dir / "capcruncher-local-apptainer" / "profile.v9+.yaml").exists()
    assert (profiles_dir / "capcruncher-slurm" / "profile.v9+.yaml").exists()
    assert (profiles_dir / "capcruncher-slurm-apptainer" / "profile.v9+.yaml").exists()
    assert not list(profiles_dir.glob("*/config.yaml"))
    assert "executor: slurm" in (
        profiles_dir / "capcruncher-slurm" / "profile.v9+.yaml"
    ).read_text()
    slurm_apptainer_profile = (
        profiles_dir / "capcruncher-slurm-apptainer" / "profile.v9+.yaml"
    ).read_text()
    assert "software-deployment-method:" in slurm_apptainer_profile
    assert "retries: 3" in slurm_apptainer_profile
    assert 'mem: "4G"' in slurm_apptainer_profile
    assert "mem_mb:" not in slurm_apptainer_profile


def test_pipeline_init_subcommand_installs_presets(cli_runner, tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))

    result = cli_runner.invoke(cli, ["pipeline", "init"])

    assert result.exit_code == 0
    assert (tmp_path / "snakemake" / "capcruncher-local" / "profile.v9+.yaml").exists()


def test_pipeline_uses_installed_preset(cli_runner, tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    init_result = cli_runner.invoke(cli, ["pipeline-init"])
    assert init_result.exit_code == 0

    recorded_calls = []

    class CompletedProcess:
        def __init__(self, returncode=0, stdout=b""):
            self.returncode = returncode
            self.stdout = stdout

    def fake_run(cmd, *args, **kwargs):
        recorded_calls.append(cmd)
        return CompletedProcess()

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = cli_runner.invoke(
        cli, ["pipeline", "run", "--preset", "capcruncher-local", "--no-logo", "-n"]
    )

    assert result.exit_code == 0
    assert len(recorded_calls) == 1
    first_call = recorded_calls[0]
    expected_profile = tmp_path / "snakemake" / "capcruncher-local"
    assert "--profile" in first_call
    assert str(expected_profile) in first_call
    assert "--cores" in first_call
    assert "1" in first_call


def test_pipeline_run_subcommand_uses_installed_preset(
    cli_runner, tmp_path, monkeypatch
):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    init_result = cli_runner.invoke(cli, ["pipeline", "init"])
    assert init_result.exit_code == 0

    recorded_calls = []

    class CompletedProcess:
        def __init__(self, returncode=0, stdout=b""):
            self.returncode = returncode
            self.stdout = stdout

    def fake_run(cmd, *args, **kwargs):
        recorded_calls.append(cmd)
        return CompletedProcess()

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = cli_runner.invoke(
        cli, ["pipeline", "run", "--preset", "capcruncher-local", "--no-logo", "-n"]
    )

    assert result.exit_code == 0
    assert len(recorded_calls) == 1
    first_call = recorded_calls[0]
    expected_profile = tmp_path / "snakemake" / "capcruncher-local"
    assert first_call[first_call.index("--profile") + 1] == str(expected_profile)
    assert "--cores" in first_call
    assert "1" in first_call


def test_pipeline_legacy_invocation_warns_with_new_command(
    cli_runner, tmp_path, monkeypatch
):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    init_result = cli_runner.invoke(cli, ["pipeline-init"])
    assert init_result.exit_code == 0

    recorded_calls = []

    class CompletedProcess:
        def __init__(self, returncode=0, stdout=b""):
            self.returncode = returncode
            self.stdout = stdout

    def fake_run(cmd, *args, **kwargs):
        recorded_calls.append(cmd)
        return CompletedProcess()

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = cli_runner.invoke(
        cli, ["pipeline", "--preset", "capcruncher-local", "--no-logo", "-n"]
    )

    assert result.exit_code == 0
    assert "Use 'capcruncher pipeline run ...' instead" in result.output
    assert recorded_calls


def test_pipeline_touches_outputs_after_real_run(cli_runner, tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    init_result = cli_runner.invoke(cli, ["pipeline-init"])
    assert init_result.exit_code == 0

    recorded_calls = []

    class CompletedProcess:
        def __init__(self, returncode=0, stdout=b""):
            self.returncode = returncode
            self.stdout = stdout

    def fake_run(cmd, *args, **kwargs):
        recorded_calls.append(cmd)
        return CompletedProcess()

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = cli_runner.invoke(
        cli, ["pipeline", "run", "--preset", "capcruncher-local", "--no-logo"]
    )

    assert result.exit_code == 0
    assert len(recorded_calls) == 2
    assert "--touch" not in recorded_calls[0]
    assert "--touch" in recorded_calls[1]


def test_pipeline_init_copies_nested_profile_files(tmp_path, monkeypatch):
    package_root = tmp_path / "package"
    source_profile = package_root / "pipeline" / "profiles" / "local"
    source_profile.mkdir(parents=True)
    (source_profile / "profile.v9+.yaml").write_text(
        "executor: local\n", encoding="utf-8"
    )
    (source_profile / "scripts").mkdir()
    (source_profile / "scripts" / "submit.sh").write_text(
        "#!/usr/bin/env bash\n", encoding="utf-8"
    )

    monkeypatch.setattr(
        cli_pipeline.resources,
        "files",
        lambda package: SimpleNamespace(joinpath=lambda *parts: package_root.joinpath(*parts)),
    )

    destination = cli_pipeline.install_pipeline_preset(
        "capcruncher-local", tmp_path / "profiles", force=False
    )

    assert (destination / "profile.v9+.yaml").exists()
    assert (destination / "scripts" / "submit.sh").exists()


def test_pipeline_accepts_legacy_preset_alias(cli_runner, tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    init_result = cli_runner.invoke(cli, ["pipeline-init"])
    assert init_result.exit_code == 0

    recorded_calls = []

    class CompletedProcess:
        def __init__(self, returncode=0, stdout=b""):
            self.returncode = returncode
            self.stdout = stdout

    def fake_run(cmd, *args, **kwargs):
        recorded_calls.append(cmd)
        return CompletedProcess()

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = cli_runner.invoke(cli, ["pipeline", "run", "--preset", "local", "--no-logo", "-n"])

    assert result.exit_code == 0
    expected_profile = tmp_path / "snakemake" / "capcruncher-local"
    assert recorded_calls[0][recorded_calls[0].index("--profile") + 1] == str(
        expected_profile
    )


def test_pipeline_preset_forwards_container_config(cli_runner, tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    init_result = cli_runner.invoke(cli, ["pipeline-init"])
    assert init_result.exit_code == 0

    recorded_calls = []

    class CompletedProcess:
        def __init__(self, returncode=0, stdout=b""):
            self.returncode = returncode
            self.stdout = stdout

    def fake_run(cmd, *args, **kwargs):
        recorded_calls.append(cmd)
        return CompletedProcess()

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = cli_runner.invoke(
        cli,
        [
            "pipeline",
            "run",
            "--preset",
            "capcruncher-local-apptainer",
            "--no-logo",
            "-n",
            "--config",
            "execution.container_image=docker://example/capcruncher:test",
        ],
    )

    assert result.exit_code == 0
    first_call = recorded_calls[0]
    expected_profile = tmp_path / "snakemake" / "capcruncher-local-apptainer"
    assert first_call[first_call.index("--profile") + 1] == str(expected_profile)
    assert "--config" in first_call
    assert "execution.container_image=docker://example/capcruncher:test" in first_call


def test_pipeline_scale_resources_sets_environment(cli_runner, tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    init_result = cli_runner.invoke(cli, ["pipeline-init"])
    assert init_result.exit_code == 0

    recorded_calls = []

    class CompletedProcess:
        def __init__(self, returncode=0, stdout=b""):
            self.returncode = returncode
            self.stdout = stdout

    def fake_run(cmd, *args, **kwargs):
        recorded_calls.append((cmd, kwargs.get("env")))
        return CompletedProcess()

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = cli_runner.invoke(
        cli,
        [
            "pipeline",
            "run",
            "--preset",
            "capcruncher-slurm-apptainer",
            "--scale-resources",
            "1.5",
            "--no-logo",
            "-n",
        ],
    )

    assert result.exit_code == 0
    first_call, first_env = recorded_calls[0]
    expected_profile = tmp_path / "snakemake" / "capcruncher-slurm-apptainer"
    assert first_call[first_call.index("--profile") + 1] == str(expected_profile)
    assert first_env["SCALE_RESOURCES"] == "1.5"


def test_pipeline_does_not_add_default_cores_for_equals_form(
    cli_runner, tmp_path, monkeypatch
):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    init_result = cli_runner.invoke(cli, ["pipeline-init"])
    assert init_result.exit_code == 0

    recorded_calls = []

    class CompletedProcess:
        def __init__(self, returncode=0, stdout=b""):
            self.returncode = returncode
            self.stdout = stdout

    def fake_run(cmd, *args, **kwargs):
        recorded_calls.append(cmd)
        return CompletedProcess()

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = cli_runner.invoke(
        cli,
        ["pipeline", "run", "--preset", "capcruncher-local", "--no-logo", "--cores=8", "-n"],
    )

    assert result.exit_code == 0
    first_call = recorded_calls[0]
    assert "--cores=8" in first_call
    assert "--cores" not in first_call


def test_pipeline_rejects_preset_and_profile_together(cli_runner, tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    init_result = cli_runner.invoke(cli, ["pipeline-init"])
    assert init_result.exit_code == 0

    result = cli_runner.invoke(
        cli,
        [
            "pipeline",
            "run",
            "--preset",
            "local",
            "--no-logo",
            "--profile=custom-profile",
            "-n",
        ],
    )

    assert result.exit_code != 0
    assert "Use either --preset or --profile" in result.output


@pytest.mark.parametrize(
    "infile,flags",
    [
        (
            "chr14.fa.gz",
            ["-r", "dpnii", "--sort"],
        ),
        pytest.param(
            "chr14.fa.gz",
            [
                "-r",
                "dpn",
            ],
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_genome_digest(cli_runner, data_pipeline, tmpdir, infile, flags):
    infile = os.path.join(data_pipeline, infile)
    outfile = os.path.join(tmpdir, "digested.bed")
    tmp_log = os.path.join(tmpdir, "gd.log")

    result = cli_runner.invoke(
        cli, ["genome", "digest", infile, "-o", outfile, "-l", tmp_log, *flags]
    )
    assert result.exit_code == 0
    assert os.path.exists(outfile)


def test_fastq_split_python(cli_runner, data_pipeline, tmpdir):
    output_prefix = pathlib.Path(tmpdir) / "split" / "sample"
    output_prefix.parent.mkdir()

    result = cli_runner.invoke(
        cli,
        [
            "fastq",
            "split",
            os.path.join(data_pipeline, "SAMPLE-A_REP1_1.fastq.gz"),
            os.path.join(data_pipeline, "SAMPLE-A_REP1_2.fastq.gz"),
            "-m",
            "python",
            "-o",
            str(output_prefix),
            "-n",
            "1000000",
            "--gzip",
        ],
    )

    assert result.exit_code == 0
    assert (output_prefix.parent / "sample_part0_1.fastq.gz").exists()
    assert (output_prefix.parent / "sample_part0_2.fastq.gz").exists()


@pytest.mark.parametrize(
    "infiles,outfile,flags",
    [
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "out",
            [],
        )
    ],
)
def test_deduplicate_fastq(
    cli_runner, data_deduplication, tmpdir, infiles, outfile, flags
):
    infiles = [os.path.join(data_deduplication, fn) for fn in infiles]
    outfile = os.path.join(tmpdir, outfile)

    result = cli_runner.invoke(
        cli,
        [
            "fastq",
            "deduplicate",
            "-1",
            infiles[0],
            "-2",
            infiles[1],
            "-o",
            outfile,
            *flags,
        ],
    )

    assert result.exit_code == 0

    fastq_deduped = glob.glob(f"{tmpdir}/*.fastq*")
    assert len(fastq_deduped) == len(infiles)


@pytest.mark.parametrize(
    "infiles,flags",
    [
        (
            ("digest_1.fastq.gz",),
            [
                "-m",
                "flashed",
                "-r",
                "dpnii",
                "--minimum_slice_length",
                "18",
            ],
        ),
        (
            ("digest_1.fastq.gz", "digest_2.fastq.gz"),
            [
                "-m",
                "pe",
                "-r",
                "dpnii",
                "--minimum_slice_length",
                "18",
            ],
        ),
        pytest.param(
            ("digest_1.fastq.gz", "digest_2.fastq.gz"),
            [
                "-m",
                "pe",
                "-r",
                "dpnii",
                "--minimum_slice_length",
                "18",
            ],
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_fastq_digest(cli_runner, data_digestion, tmpdir, infiles, flags):
    infiles = [os.path.join(data_digestion, fn) for fn in infiles]
    outfile = os.path.join(tmpdir, "digested.fastq")

    result = cli_runner.invoke(
        cli, ["fastq", "digest", *infiles, "-o", outfile, *flags]
    )
    assert result.exit_code == 0
    assert os.path.exists(outfile)


@pytest.mark.parametrize(
    "bam,beds,flags",
    [
        (
            "test.pe.bam",
            [
                "test_capture.bed",
            ],
            [
                "-n",
                "TEST_OVERLAP",
                "-a",
                "get",
                "-f",
                1e-9,
            ],
        ),
        (
            "test.pe.bam",
            [
                "test_capture.bed",
            ],
            [
                "-n",
                "TEST_OVERLAP",
                "-a",
                "get",
                "-f",
                1e-9,
                "--priority-chroms",
                "chr14",
                "--prioritize-cis-slices",
            ],
        ),
    ],
)
def test_alignment_annotation(cli_runner, data_annotation, tmpdir, bam, beds, flags):
    bam = os.path.join(data_annotation, bam)
    beds = [os.path.join(data_annotation, bed) for bed in beds]
    blacklist = os.path.join(data_annotation, "test_exlcusions_corrected.bed")
    outfile = os.path.join(tmpdir, "annotated.parquet")

    result = cli_runner.invoke(
        cli,
        [
            "alignments",
            "annotate",
            bam,
            "-b",
            *beds,
            "-o",
            outfile,
            *flags,
            "--blacklist",
            blacklist,
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(outfile)


@pytest.mark.parametrize(
    "mode,bam,annotations,flags",
    [
        (
            "capture",
            "test.flashed.bam",
            "test.annotations.parquet",
            [],
        ),
        (
            "tri",
            "test.flashed.bam",
            "test.annotations.parquet",
            [],
        ),
        (
            "tiled",
            "test.flashed.bam",
            "test.annotations.parquet",
            [],
        ),
    ],
)
def test_alignment_filter(
    cli_runner, data_filter, tmpdir, mode, bam, annotations, flags
):
    bam = os.path.join(data_filter, bam)
    annotations = os.path.join(data_filter, annotations)
    output_prefix = os.path.join(tmpdir, "filtered.parquet")
    stats_prefix = os.path.join(tmpdir, "test_filtering_cli")

    result = cli_runner.invoke(
        cli,
        [
            "alignments",
            "filter",
            mode,
            "-b",
            bam,
            "-a",
            annotations,
            "-o",
            output_prefix,
            "--sample-name",
            "test",
            *flags,
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0
    assert os.path.exists(f"{output_prefix}.slices.parquet")


@pytest.mark.parametrize(
    "slices,read_type,output",
    [
        (
            "slices.flashed.parquet",
            "flashed",
            "deduplicated.flashed.parquet",
        ),
        (
            "slices.pe.parquet",
            "pe",
            "deduplicated.pe.parquet",
        ),
    ],
)
def test_interactions_deduplicate(
    cli_runner, data_deduplication_alignments, tmpdir, slices, read_type, output
):
    slices = os.path.join(data_deduplication_alignments, slices)
    output = os.path.join(tmpdir, output)

    """
    slices: os.PathLike,
    output: os.PathLike,
    read_type: str = "flashed",
    sample_name: str = "sampleX",
    statistics: os.PathLike = "deduplication_stats.json",
    """

    result = cli_runner.invoke(
        cli,
        [
            "interactions",
            "deduplicate",
            slices,
            "--read-type",
            read_type,
            "-o",
            output,
            "--sample-name",
            "test",
            "--statistics",
            "test_deduplication_cli.json",
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(output)


@pytest.mark.parametrize(
    "slices,bins,viewpoints,output,flags",
    [
        ("slices.flashed.parquet", "bins.bed.gz", "viewpoints.bed", "counts.hdf5", []),
    ],
)
def test_reporters_count(
    cli_runner,
    data_deduplication_alignments,
    data_reporters_count,
    data_pipeline,
    tmpdir,
    slices,
    bins,
    viewpoints,
    output,
    flags,
):
    slices = os.path.join(data_deduplication_alignments, slices)
    bins = os.path.join(data_reporters_count, bins)
    viewpoints = os.path.join(data_pipeline, "mm9_capture_viewpoints_Slc25A37.bed")
    output = os.path.join(tmpdir, output)

    result = cli_runner.invoke(
        cli,
        [
            "interactions",
            "count",
            slices,
            "-o",
            output,
            "-f",
            bins,
            "-v",
            viewpoints,
            *flags,
        ],
    )

    logger.debug(result.stdout)
    logger.debug(result.exception)
    assert result.exit_code == 0
    assert os.path.exists(output)


def test_reporters_count_fixture_matches_viewpoint_file(
    data_reporters_count, data_pipeline
):
    reporters = os.path.join(
        os.path.dirname(data_reporters_count),
        "reporter_count",
        "SAMPLE-A_REP1.parquet",
    )
    viewpoints = os.path.join(data_pipeline, "mm9_capture_viewpoints_Slc25A37.bed")
    viewpoint_names = pd.read_csv(viewpoints, sep="\t", header=None)[3].to_list()

    assert viewpoint_names == ["Slc25A37"]

    for parquet_file in sorted(pathlib.Path(reporters).glob("*.parquet")):
        reporter_viewpoints = pd.read_parquet(parquet_file, columns=["viewpoint"])[
            "viewpoint"
        ]
        assert reporter_viewpoints.cat.categories.to_list() == viewpoint_names
        assert reporter_viewpoints.value_counts(dropna=False).to_dict() == {
            "Slc25A37": len(reporter_viewpoints)
        }


@pytest.mark.parametrize(
    "command,cooler_fn,bin_size,output,flags",
    [
        ("bin", "SAMPLE-A_REP1.hdf5", int(1e5), "binned.hdf5", []),
        (
            "fragments-to-bins",
            "SAMPLE-A_REP1.hdf5",
            int(1e5),
            "binned_legacy.hdf5",
            [],
        ),
    ],
)
def test_reporters_store_binned(
    cli_runner,
    data_reporters_store,
    tmpdir,
    command,
    cooler_fn,
    bin_size,
    output,
    flags,
):
    clr = os.path.join(data_reporters_store, cooler_fn)
    output = os.path.join(tmpdir, output)

    result = cli_runner.invoke(
        cli,
        [
            "interactions",
            command,
            clr,
            "-o",
            output,
            "-b",
            bin_size,
            *flags,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(output)


def test_countable_reporters_only_include_bed_viewpoint_categories(tmp_path):
    viewpoints = tmp_path / "viewpoints.bed"
    reporters = tmp_path / "reporters.parquet"
    output = tmp_path / "countable"

    viewpoints.write_text("chr14\t69902454\t69903469\tSlc25A37\n")
    pd.DataFrame(
        {
            "viewpoint": pd.Categorical(
                ["Slc25A37", "Slc25A37", "Slc25A37"],
                categories=["Slc25A37", "reporters_pe_80", "duplicate_coords_1"],
            ),
            "parent_id": [1, 2, 3],
            "restriction_fragment": [10, 20, 30],
        }
    ).to_parquet(reporters)

    cleaned = _write_countable_reporters(reporters, viewpoints, output)
    cleaned_df = pd.read_parquet(cleaned)

    assert cleaned_df["viewpoint"].cat.categories.to_list() == ["Slc25A37"]
    assert cleaned_df["viewpoint"].to_list() == ["Slc25A37"] * 3


def test_countable_reporters_reject_actual_non_viewpoint_values(tmp_path):
    viewpoints = tmp_path / "viewpoints.bed"
    reporters = tmp_path / "reporters.parquet"
    output = tmp_path / "countable"

    viewpoints.write_text("chr14\t69902454\t69903469\tSlc25A37\n")
    pd.DataFrame(
        {
            "viewpoint": ["Slc25A37", "reporters_pe_80"],
            "parent_id": [1, 2],
            "restriction_fragment": [10, 20],
        }
    ).to_parquet(reporters)

    with pytest.raises(ValueError, match="reporters_pe_80"):
        _write_countable_reporters(reporters, viewpoints, output)


@pytest.mark.parametrize(
    "infiles,viewpoint,output,flags",
    [
        (["SAMPLE-A_REP1.hdf5", "SAMPLE-A_REP1.hdf5"], "Slc25A37", "merged.hdf5", []),
    ],
)
def test_reporters_store_merge(
    cli_runner,
    data_reporters_store,
    tmpdir,
    infiles,
    viewpoint,
    output,
    flags,
):
    import cooler

    hdf5_files = [os.path.join(data_reporters_store, fn) for fn in infiles]
    output = os.path.join(tmpdir, output)

    result = cli_runner.invoke(
        cli,
        [
            "interactions",
            "merge",
            *hdf5_files,
            "-o",
            output,
            *flags,
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(output)
    assert (
        cooler.Cooler(f"{output}::{viewpoint}").pixels()[:]["count"].sum()
        == cooler.Cooler(f"{hdf5_files[0]}::{viewpoint}").pixels()[:]["count"].sum() * 2
    )


@pytest.mark.parametrize(
    "cooler_fn,output_prefix,outfile,flags",
    [
        ("SAMPLE-A_REP1.hdf5", "test", "test_Slc25A37.bedgraph", []),
        ("SAMPLE-A_REP1.hdf5", "test", "test_Slc25A37.bigWig", ["-f", "bigwig"]),
    ],
)
def test_reporters_pileup(
    cli_runner,
    data_reporters_store,
    tmpdir,
    cooler_fn,
    output_prefix,
    outfile,
    flags,
):
    clr = os.path.join(data_reporters_store, cooler_fn)
    output = os.path.join(tmpdir, output_prefix)
    outfile = os.path.join(tmpdir, outfile)

    result = cli_runner.invoke(
        cli,
        [
            "interactions",
            "pileup",
            clr,
            "-o",
            output,
            *flags,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(outfile)


@pytest.fixture
def fragments_file(tmpdir):
    fragments = os.path.join(tmpdir, "fragments.bed")
    with open(fragments, "w") as file:
        file.write("chr1\t100\t200\tfragment1\n")
        file.write("chr2\t300\t400\tfragment2\n")
    return fragments


@pytest.fixture
def viewpoints_file(tmpdir):
    viewpoints = os.path.join(tmpdir, "viewpoints.bed")
    with open(viewpoints, "w") as file:
        file.write("chr1\t150\t160\tviewpoint1\n")
        file.write("chr2\t350\t360\tviewpoint2\n")
    return viewpoints


def test_make_chicago_maps(cli_runner, tmpdir, fragments_file, viewpoints_file):
    outputdir = str(tmpdir)

    result = cli_runner.invoke(
        cli,
        [
            "utilities",
            "make-chicago-maps",
            "--fragments",
            fragments_file,
            "--viewpoints",
            viewpoints_file,
            "-o",
            outputdir,
        ],
    )

    assert result.exit_code == 0

    # Check if the renamed fragments file exists
    fragments_new = os.path.join(outputdir, "fragments.rmap")
    assert os.path.exists(fragments_new)

    # Check if the baitmap file exists and has the correct content
    baitmap_file = os.path.join(outputdir, "viewpoints.baitmap")
    assert os.path.exists(baitmap_file)

    with open(baitmap_file, "r") as file:
        content = file.read()
        assert "chr1\t100\t200\tfragment1\tviewpoint" in content


def test_cis_and_trans_stats_accepts_empty_parquet_directory(cli_runner, tmp_path):
    slices = tmp_path / "empty_slices.parquet"
    slices.mkdir()
    output = tmp_path / "cis_and_trans.json"

    result = cli_runner.invoke(
        cli,
        [
            "utilities",
            "cis-and-trans-stats",
            str(slices),
            "--assay",
            "tiled",
            "--sample-name",
            "SAMPLE-A",
            "-o",
            str(output),
        ],
    )

    assert result.exit_code == 0
    assert output.read_text() == '{"stats":[]}'
