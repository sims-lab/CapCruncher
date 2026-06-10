# CapCruncher Agent Notes

## Repository Workflow

- Work on the local `develop` branch unless the user explicitly asks otherwise.
- Use Conventional Commits, for example `refactor: migrate intervals to pyranges1`.
- Preserve unrelated dirty workspace files. Recent local-only files such as `plans/`
  and `.gitignore` entries may be user context rather than repo changes.
- Prefer narrow, tested changes. This repo has been unmaintained for a while, so
  fragile modernization fixes should usually get regression tests.

## Runtime Targets

- Target Python 3.12+ only.
- The supported validation environment is conda env `cc`:
  - `conda run -n cc pytest ...`
  - `conda run -n cc python -m py_compile ...`
- The repo `.venv` is Python 3.12 but may not include `pytest`.
- Treat `capcruncher-tools` as an external dependency. Keep
  `capcruncher/cli/interactions_count.py` delegating to it unless the dependency
  is proven broken after updating.
- `capcruncher-tools` should be at the latest Polars-compatible release in
  requirements and workflow envs.
- Reporter counting should decide whether to call `capcruncher-tools` from a
  quick row-count summary. Viewpoint categories with zero reporter rows should
  warn and skip cooler creation rather than writing empty cooler groups.

## Modernization Constraints

- Use `pyranges1` only: `import pyranges1 as pr`.
- Do not add fallback support for original `pyranges`.
- `pyranges1.PyRanges` is a DataFrame subclass. Use DataFrame operations
  directly instead of old `.df` or `.as_df()` accessors.
- `pybedtools` is removed from active code and dependency manifests. Use
  `pyranges1`, pandas, polars, or pysam instead.
- `ibis` has been removed. Use polars or direct DuckDB/pandas as appropriate.
- Assume active CapCruncher code remains PyRanges1-only.

## Pipeline And Snakemake

- Snakemake 9 is the target.
- Pipeline presets are installed to the standard Snakemake user profile path:
  `${XDG_CONFIG_HOME:-~/.config}/snakemake`.
- Bundled preset names are namespaced:
  - `capcruncher-local`
  - `capcruncher-local-conda`
  - `capcruncher-local-apptainer`
  - `capcruncher-slurm`
  - `capcruncher-slurm-apptainer`
- Legacy short aliases such as `local` and `slurm-apptainer` are accepted for
  compatibility, but docs should prefer the namespaced preset names.
- `pipeline-init` should copy profile directories recursively. Future profile
  helper scripts or nested files must not be silently dropped.
- `capcruncher pipeline -n` / `--dry-run` must not trigger the follow-up
  `snakemake --touch`.
- `--scale-resources` is CapCruncher-specific. It sets `SCALE_RESOURCES` for
  workflow resource functions; it is not a native Snakemake setting.
- Pipeline tests are marked with `pipeline` and `slow`. CI runs quick tests with
  `-m "not pipeline"` and then runs the slow pipeline suite separately with
  `-m "pipeline"`.
- When running pipeline tests, pass `--cores 4` unless there is a specific reason
  to reproduce single-core behavior.

## Containers

- Docker image builds should support `linux/amd64` and `linux/arm64`.
- macOS users run Linux containers through Docker Desktop or Colima. Apple
  Silicon should use `linux/arm64` automatically, with `--platform` available
  when needed.
- The Docker image includes Apptainer so it can run containerised Snakemake
  workflows from inside Docker where host privileges allow it.
- Apptainer is the supported runtime on HPC. Docker is for local workstation and
  CI usage, not shared HPC execution.
- Container smoke checks should include:
  - `capcruncher --help`
  - `apptainer --version`
  - `quarto --version`

## Plotting

- Plotting should use PlotNado, not the removed CapCruncher CoolBox plotting API.
- The pipeline writes PlotNado TOML templates alongside generated figures.
- Advanced plotting docs/notebooks should show this workflow:
  1. Run the pipeline with `plot.create: True`.
  2. Locate `capcruncher_output/results/plots/templates/*.toml`.
  3. Edit the TOML or load it with `plotnado.GenomicFigure.from_toml`.
  4. Render with `capcruncher plot --template ... --region ... --output ...`
     or `fig.save(...)`.
- Keep examples clean and runnable with small test data where possible.

## UCSC Hub / TrackNado

- `make_ucsc_hub.py` uses TrackNado `HubBuilder` with a metadata extractor.
- Track metadata parsing should be strict and fail clearly for unsupported
  CapCruncher output filenames.
- Custom UCSC hub genomes require a twoBit file. Fail early with a clear error
  if `custom_genome` is true and `genome.twobit` is missing.
- Do not keep unused hub scaffold keys unless TrackNado actually consumes them.

## Useful Checks

Scan for legacy interval/query dependencies:

```bash
rg -n "import pyranges as pr|from pyranges\\b|pybedtools|BedTool|pyranges<=|\\.as_df\\(|pyranges_to_dataframe|ibis-framework|\\bibis\\b" \
  capcruncher tests requirements*.txt environment.yml capcruncher/pipeline/workflow/envs/environment.yml
```

Focused modernization checks:

```bash
conda run -n cc pytest tests/test_cli.py -q
conda run -n cc pytest tests/test_workflow_scripts.py tests/test_plotting.py -q
conda run -n cc pytest tests/test_annotate.py tests/test_pileup.py tests/test_storage_api.py -q
conda run -n cc pytest tests/test_utils.py -q -k 'bed_validation_and_formatting or intersect_bins_and_interval_conversion or read_dataframes'
conda run -n cc pytest -q -m "not pipeline"
conda run -n cc pytest -q -m "pipeline" --cores 4
conda run -n cc python -m py_compile capcruncher/api/storage.py capcruncher/cli/cli_pipeline.py capcruncher/pipeline/workflow/scripts/make_ucsc_hub.py capcruncher/utils.py
```

Known environment caveat:

- `tests/test_utils.py::test_viewpoint_coordinates` can be blocked by
  `ModuleNotFoundError: capcruncher_tools.capcruncher_tools`; treat that as an
  environment/dependency issue unless changing that path.

## Recent Relevant Commits

- `c822b22 chore: add plans directory to .gitignore`
- `23987f0 test: cover brittle pipeline modernization paths`
- `9a46b17 feat: modernize pipeline profiles and container docs`
- `3c224ba docs: point plotting customization to plotnado`
- `213f5cf refactor: align hub generation with tracknado extractors`
- `96a3381 chore: prune unused dependencies`
