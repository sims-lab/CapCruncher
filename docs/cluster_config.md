# Set Up Snakemake Execution Presets

CapCruncher ships Snakemake 9 execution presets that use executor plugins rather than legacy `cluster`, `drmaa`, or submit-script profiles.

Apptainer is the supported container runtime for HPC execution. Docker is useful
for local and CI runs, but shared HPC systems generally require Apptainer.

Install the bundled presets with:

```bash
capcruncher pipeline init
```

By default this writes profiles to:

```text
${XDG_CONFIG_HOME:-~/.config}/snakemake
```

This is Snakemake's standard user profile directory on Linux, so users can edit
the generated profile files directly and Snakemake can also find them by name.

## Edit Profiles

Each preset is a normal Snakemake profile directory containing `profile.v9+.yaml`.
For example, edit the local Apptainer profile with:

```bash
${EDITOR:-nano} "${XDG_CONFIG_HOME:-$HOME/.config}/snakemake/capcruncher-local-apptainer/profile.v9+.yaml"
```

Common edits include `jobs`, `latency-wait`, `retries`, `default-resources`,
`apptainer-args`, and SLURM executor options. Values in the profile become
Snakemake defaults. Command-line options passed to `capcruncher pipeline run` still
take precedence for that run.

To create a site-specific profile without losing future bundled defaults, copy
one of the installed profiles and edit the copy:

```bash
cp -r \
  "${XDG_CONFIG_HOME:-$HOME/.config}/snakemake/capcruncher-slurm-apptainer" \
  "${XDG_CONFIG_HOME:-$HOME/.config}/snakemake/my-lab-capcruncher"

${EDITOR:-nano} "${XDG_CONFIG_HOME:-$HOME/.config}/snakemake/my-lab-capcruncher/profile.v9+.yaml"
capcruncher pipeline run --preset my-lab-capcruncher --jobs 50
```

## Update Profiles

Running `capcruncher pipeline init` again will not overwrite existing profiles.
This protects local edits. To refresh the installed CapCruncher defaults after
upgrading CapCruncher, use `--force`:

```bash
capcruncher pipeline init --force
```

To refresh only one bundled preset:

```bash
capcruncher pipeline init --preset capcruncher-slurm-apptainer --force
```

`--force` replaces the selected installed profile directories. Back up any local
edits first, or keep site-specific changes in a copied profile such as
`my-lab-capcruncher`.

Run the pipeline with a preset:

```bash
capcruncher pipeline run --preset capcruncher-local -n
capcruncher pipeline run --preset capcruncher-slurm --jobs 50
capcruncher pipeline run --preset capcruncher-slurm-apptainer --jobs 50
```

The bundled SLURM preset is a Snakemake 9 profile:

```yaml
executor: slurm
jobs: 100
latency-wait: 60
printshellcmds: true
rerun-incomplete: true
retries: 3
show-failed-logs: true
default-resources:
  mem: "4G"
  runtime: 60
```

For Apptainer execution, use the `capcruncher-slurm-apptainer` preset:

```yaml
executor: slurm
jobs: 100
latency-wait: 60
software-deployment-method:
  - apptainer
apptainer-args: --cleanenv
printshellcmds: true
rerun-incomplete: true
retries: 3
show-failed-logs: true
default-resources:
  mem: "4G"
  runtime: 60
```

Modern Snakemake separates SLURM executor options from SLURM-specific
resources. Executor options control submission behavior, status checks, and log
handling. For example, QoS, reservation, and SLURM log directory can be passed
through the SLURM executor plugin:

```bash
capcruncher pipeline run \
  --preset capcruncher-slurm \
  --slurm-qos normal \
  --slurm-reservation reservation-name \
  --slurm-logdir logs/slurm \
  --jobs 50
```

SLURM account and partition are resources named `slurm_account` and
`slurm_partition`. Set them with `--default-resources`, rule-specific
`--set-resources`, or directly in the editable profile:

```bash
capcruncher pipeline run \
  --preset capcruncher-slurm \
  --default-resources \
      slurm_account=my-account \
      slurm_partition=standard \
      mem=6G \
      runtime=120 \
      disk_mb=100000 \
  --set-resources \
      align_bowtie2:mem=12G \
      align_bowtie2:runtime=240 \
      split:slurm_partition=long
```

For persistent cluster defaults, edit the installed profile instead of repeating
these values on every command:

```yaml
executor: slurm
jobs: 100
slurm-logdir: logs/slurm
slurm-qos: normal
default-resources:
  mem: "6G"
  runtime: 120
  disk_mb: 100000
  slurm_account: "my-account"
  slurm_partition: "standard"
set-resources:
  align_bowtie2:
    mem: "12G"
    runtime: 240
  split:
    slurm_partition: "long"
```

The SLURM executor plugin can also select partitions automatically from a
partition limits file via `--slurm-partition-config`. That is usually better
than hard-coding partitions when a cluster has multiple queues with different
runtime, memory, or CPU limits.

CapCruncher's `--scale-resources` is a convenience wrapper around the workflow's
dynamic `mem` and `runtime` resources. It multiplies CapCruncher-authored memory
and runtime requests, while Snakemake still applies retry-aware resource
functions via the rule `attempt`. Use it when most jobs need a global uplift
without editing the profile:

```bash
capcruncher pipeline run \
  --preset capcruncher-slurm-apptainer \
  --scale-resources 1.5 \
  --jobs 50
```
