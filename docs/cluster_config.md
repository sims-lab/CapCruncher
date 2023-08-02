
# Set-up a Snakemake profile

This is not essential but it will make running the pipeline much easier by submitting jobs to the cluster automatically and using pre-set parameters.

**Note:** Cookiecutter is required for this step. This can be installed using `pip install cookiecutter`.


### For SLURM based clusters:

```bash
# create config directory that snakemake searches for profiles (or use something else)
profile_dir="${HOME}/.config/snakemake"
mkdir -p "$profile_dir"
# use cookiecutter to create the profile in the config directory
template="gh:Snakemake-Profiles/slurm"
cookiecutter --output-dir "$profile_dir" "$template"
```

### For SGE based clusters:

!!! warning
    This has not been tested

```bash
mkdir -p ~/.config/snakemake
cd ~/.config/snakemake
cookiecutter https://github.com/Snakemake-Profiles/sge.git
```

### Example SLURM profile:

```
/home/a/asmith/.config/snakemake/slurm/
├── config.yaml
├── CookieCutter.py
├── __pycache__
│   ├── CookieCutter.cpython-310.pyc
│   ├── CookieCutter.cpython-311.pyc
│   ├── slurm_utils.cpython-310.pyc
│   └── slurm_utils.cpython-311.pyc
├── settings.json
├── slurm-jobscript.sh
├── slurm-sidecar.py
├── slurm-status.py
├── slurm-submit.py
└── slurm_utils.py
```

`settings.json`:

```json
{
    "SBATCH_DEFAULTS": "--partition=short --time=0-01:00:00 --mem=3G",
    "CLUSTER_NAME": "",
    "CLUSTER_CONFIG": ""
}
```

`config.yaml`:

```yaml

cluster-sidecar: "slurm-sidecar.py"
cluster-cancel: "scancel"
restart-times: "0"
jobscript: "slurm-jobscript.sh"
cluster: "slurm-submit.py"
cluster-status: "slurm-status.py"
max-jobs-per-second: "10"
max-status-checks-per-second: "10"
local-cores: 1
latency-wait: "5"
use-conda: "True"
use-singularity: "False"
singularity-args: -B /ceph  -B /databank -B $TMPDIR --cleanenv
jobs: "50"
printshellcmds: "True"
retries: 3

# Example resource configuration
# default-resources:
#   - runtime=100
#   - mem_mb=6000
#   - disk_mb=1000000
# # set-threads: map rule names to threads
# set-threads:
#   - single_core_rule=1
#   - multi_core_rule=10
# # set-resources: map rule names to resources in general
# set-resources:
#   - high_memory_rule:mem_mb=12000
#   - long_running_rule:runtime=1200

```

**Note**: The singularity-args are required to mount the data directories into the container. e.g.

```bash
singularity-args: -B /ceph  -B /databank
```

Gives the container access to the `/ceph` and `/databank` directories on the cluster. The current working directory is also mounted into the container by default. You can add additional directories by adding more `-B` flags. Obviously this will be different for each cluster so you'll need your own defaults. The `$TMPDIR` is also mounted as this causes errors if not. The `--cleanenv` flag is also required to prevent the container from inheriting the environment from the host.
