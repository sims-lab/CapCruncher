# syntax=docker/dockerfile:1

# ── Stage 1: install conda env + build tools ──────────────────────────────────
# Separated so environment.yml changes don't bust the source-code layer and
# vice-versa.
FROM mambaorg/micromamba:2.0.5 AS builder

ARG CAPCRUNCHER_VERSION=0.0.0+container
ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/capcruncher-environment.yml

# Cache conda package tarballs across builds; clean flag not needed because
# pkgs dir is a BuildKit cache mount (never written to the layer).
RUN --mount=type=cache,target=/opt/conda/pkgs,uid=1000,gid=1000,sharing=locked \
    micromamba install -y -n base -c conda-forge -c bioconda \
        apptainer \
        cxx-compiler \
        rust && \
    micromamba install -y -n base -f /tmp/capcruncher-environment.yml

WORKDIR /opt/capcruncher
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md LICENSE MANIFEST.in ./
COPY --chown=$MAMBA_USER:$MAMBA_USER capcruncher ./capcruncher

RUN --mount=type=cache,target=/home/mambauser/.cache/pip,uid=1000,gid=1000,sharing=locked \
    ln -sf /opt/conda/bin/flash2 /opt/conda/bin/flash && \
    printf 'setuptools<80\n' > /tmp/pip-build-constraints.txt && \
    PIP_CACHE_DIR=/home/mambauser/.cache/pip \
    PIP_CONSTRAINT=/tmp/pip-build-constraints.txt \
    SETUPTOOLS_SCM_PRETEND_VERSION_FOR_CAPCRUNCHER="${CAPCRUNCHER_VERSION}" \
        pip install --no-deps . && \
    micromamba remove -y -n base c-compiler cxx-compiler gcc gxx rust && \
    micromamba clean --all --yes && \
    find /opt/conda -type f \( -name "*.pyc" -o -name "*.pyo" -o -name "*.a" \) -delete && \
    find /opt/conda -type d \( -name "__pycache__" -o -name "tests" -o -name "test" \) \
        -prune -exec rm -rf '{}' + && \
    find /opt/capcruncher -type d \( -name "__pycache__" -o -name "*.egg-info" \) \
        -prune -exec rm -rf '{}' + && \
    rm -rf \
        /opt/capcruncher/dist \
        /tmp/pip-build-constraints.txt

# ── Stage 2: runtime image ────────────────────────────────────────────────────
FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.title="CapCruncher" \
      org.opencontainers.image.source="https://github.com/sims-lab/CapCruncher" \
      org.opencontainers.image.licenses="GPL-3.0-only"

COPY --link --from=builder /opt/conda      /opt/conda
COPY --link --from=builder /opt/capcruncher /opt/capcruncher

ENV CONDA_PREFIX=/opt/conda \
    MPLCONFIGDIR=/tmp/matplotlib \
    PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    PYTHONUNBUFFERED=1 \
    XDG_CACHE_HOME=/tmp/.cache

WORKDIR /opt/capcruncher

ENTRYPOINT ["capcruncher"]
CMD ["--help"]
