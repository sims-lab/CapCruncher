# Docker Usage

CapCruncher publishes a Docker/OCI image containing the CLI, Snakemake workflow, native command-line tools, Python dependencies, Apptainer, and the Quarto report runtime.

```bash
docker pull ghcr.io/sims-lab/capcruncher:latest
docker run --rm ghcr.io/sims-lab/capcruncher:latest --help
```

## macOS

Docker runs Linux containers on macOS through Docker Desktop or an equivalent
runtime such as Colima. The published CapCruncher image is multi-platform for
`linux/amd64` and `linux/arm64`, so Intel Macs and Apple Silicon Macs can use
the same image name:

```bash
docker run --rm ghcr.io/sims-lab/capcruncher:latest --help
```

On Apple Silicon, Docker should select the `linux/arm64` image automatically.
To force a specific platform, pass `--platform`:

```bash
docker run --rm --platform linux/arm64 ghcr.io/sims-lab/capcruncher:latest --help
docker run --rm --platform linux/amd64 ghcr.io/sims-lab/capcruncher:latest --help
```

When mounting project or reference data on macOS, make sure the host paths are
shared with Docker Desktop or Colima before running the pipeline.

## Run the Pipeline in Docker

Run Docker from the project directory containing `capcruncher_config.yml` and the FASTQ files or symlinks that the pipeline should process:

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v "$PWD":/work \
  -w /work \
  ghcr.io/sims-lab/capcruncher:latest \
  pipeline --cores 8
```

This runs the entire `capcruncher pipeline` command inside the Docker container and writes `capcruncher_output/` back into the mounted working directory.

If your input files are symlinks to paths outside the current directory, mount those external paths too. Docker cannot read host paths that are not mounted into the container.

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v "$PWD":/work \
  -v /data/reference:/data/reference:ro \
  -v /data/fastqs:/data/fastqs:ro \
  -w /work \
  ghcr.io/sims-lab/capcruncher:latest \
  pipeline --cores 8
```

## Run Other CLI Commands

The image entrypoint is `capcruncher`, so pass CapCruncher subcommands directly after the image name:

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v "$PWD":/work \
  -w /work \
  ghcr.io/sims-lab/capcruncher:latest \
  pipeline-config
```

To open a shell inside the image:

```bash
docker run --rm -it \
  --entrypoint bash \
  -v "$PWD":/work \
  -w /work \
  ghcr.io/sims-lab/capcruncher:latest
```

## Build the Image Locally

For development or local validation:

```bash
docker build -t capcruncher:dev .
docker run --rm capcruncher:dev --help
docker run --rm --entrypoint quarto capcruncher:dev --version
```

Build a specific platform when testing macOS compatibility:

```bash
docker buildx build --platform linux/arm64 -t capcruncher:dev-arm64 --load .
docker buildx build --platform linux/amd64 -t capcruncher:dev-amd64 --load .
```

## Containerised Workflows Inside Docker

The Docker image includes `apptainer`, so it can invoke Snakemake profiles such as `capcruncher-local-apptainer` from inside the container. This is useful when you want the outer Docker image to provide CapCruncher and Snakemake, while Snakemake still executes workflow jobs with the same `docker://ghcr.io/sims-lab/capcruncher:latest` container URI used on HPC systems.

Nested container execution depends on the host Docker configuration. On many systems, Apptainer inside Docker requires additional privileges and cache mounts:

```bash
docker run --rm -it \
  --privileged \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -e APPTAINER_CACHEDIR=/work/.apptainer-cache \
  -v "$PWD":/work \
  -w /work \
  ghcr.io/sims-lab/capcruncher:latest \
  pipeline --preset capcruncher-local-apptainer --cores 8
```

If nested Apptainer is not available on your Docker host, run `pipeline --cores 8` inside Docker instead, or run CapCruncher directly on the host with the `capcruncher-local-apptainer` or `capcruncher-slurm-apptainer` presets.

## Apptainer on HPC

Apptainer is the supported container runtime on HPC systems. Docker is intended for local workstation and CI usage; shared clusters generally do not expose a Docker daemon to users.

You can run the CapCruncher image directly with Apptainer:

```bash
apptainer exec docker://ghcr.io/sims-lab/capcruncher:latest capcruncher --help
```

From a project directory containing `capcruncher_config.yml` and input FASTQs:

```bash
apptainer exec \
  --bind "$PWD":/work \
  --pwd /work \
  docker://ghcr.io/sims-lab/capcruncher:latest \
  capcruncher pipeline --cores 8
```

For cluster-scale runs, install the editable Snakemake profiles and use the Apptainer preset so each workflow job runs through Snakemake's supported container backend:

Use the `capcruncher-local-apptainer` or `capcruncher-slurm-apptainer` presets when you want Snakemake to execute workflow jobs through its supported container deployment backend. Those presets use the same container image via the `docker://ghcr.io/sims-lab/capcruncher:latest` URI configured in `capcruncher_config.yml`.

```bash
capcruncher pipeline --preset capcruncher-local-apptainer --cores 8
capcruncher pipeline --preset capcruncher-slurm-apptainer --jobs 50
```

## Docker vs Apptainer

Use Docker when you want to run the whole CapCruncher CLI in one container on a workstation or CI runner.

Use Apptainer on HPC systems, either by running the image directly with `apptainer exec` or by using the `capcruncher-local-apptainer` and `capcruncher-slurm-apptainer` presets.
