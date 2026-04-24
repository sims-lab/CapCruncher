# Set Up Snakemake Execution Presets

CapCruncher ships Snakemake 9 execution presets that use executor plugins rather than legacy `cluster`, `drmaa`, or submit-script profiles.

Install the bundled presets with:

```bash
capcruncher pipeline-init
```

By default this writes profiles to:

```text
${XDG_CONFIG_HOME:-~/.config}/capcruncher/profiles
```

Run the pipeline with a preset:

```bash
capcruncher pipeline --preset local -n
capcruncher pipeline --preset slurm --jobs 50
capcruncher pipeline --preset slurm-apptainer --jobs 50
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

For Apptainer execution, use the `slurm-apptainer` preset:

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

Cluster-specific settings such as QoS, reservation, account, partition, and log directory should be passed through the SLURM executor plugin options, for example:

```bash
capcruncher pipeline \
  --preset slurm \
  --slurm-qos normal \
  --slurm-logdir logs/slurm \
  --jobs 50
```

Use `--set-resources` or `--default-resources` for rule-level tuning:

```bash
capcruncher pipeline \
  --preset slurm \
  --default-resources mem=6G runtime=120 disk_mb=100000 \
  --set-resources align_bowtie2:mem=12G split:runtime=240
```

CapCruncher also supports SeqNado-style resource scaling for retry-prone
runs. The bundled presets set `retries: 3`; on each retry, workflow resources
that use `mem` and `runtime` are recalculated from the Snakemake `attempt`
number. Use `--scale-resources` to raise all workflow-defined memory and runtime
requests without editing the profile:

```bash
capcruncher pipeline --preset slurm-apptainer --scale-resources 1.5 --jobs 50
```
