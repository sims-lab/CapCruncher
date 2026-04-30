FROM mambaorg/micromamba:2.0.5

ARG CAPCRUNCHER_VERSION=0.0.0+container
ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/capcruncher-environment.yml

WORKDIR /opt/capcruncher

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/capcruncher

ENV CARGO_HOME=/tmp/cargo \
    PIP_CACHE_DIR=/tmp/pip-cache \
    PIP_NO_CACHE_DIR=1 \
    RUSTUP_HOME=/tmp/rustup \
    XDG_CACHE_HOME=/tmp/.cache

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        apptainer \
        cxx-compiler \
        rust && \
    micromamba install -y -n base -f /tmp/capcruncher-environment.yml && \
    ln -sf /opt/conda/bin/flash2 /opt/conda/bin/flash && \
    apptainer --version && \
    printf 'setuptools<80\n' > /tmp/pip-build-constraints.txt && \
    python -m pip install --no-cache-dir --upgrade pip && \
    PIP_CONSTRAINT=/tmp/pip-build-constraints.txt \
    SETUPTOOLS_SCM_PRETEND_VERSION_FOR_CAPCRUNCHER="${CAPCRUNCHER_VERSION}" \
        python -m pip install --no-cache-dir --no-deps . && \
    micromamba remove -y -n base \
        c-compiler \
        cxx-compiler \
        gcc \
        gxx \
        rust && \
    micromamba clean --all --yes && \
    find /opt/conda -type d \( -name "__pycache__" -o -name "tests" -o -name "test" \) -prune -exec rm -rf '{}' + && \
    find /opt/conda -type f \( -name "*.pyc" -o -name "*.pyo" -o -name "*.a" \) -delete && \
    find /opt/capcruncher -type d \( -name "__pycache__" -o -name "*.egg-info" \) -prune -exec rm -rf '{}' + && \
    find /opt/capcruncher -type f \( -name "*.pyc" -o -name "*.pyo" \) -delete && \
    rm -rf \
        /opt/conda/.cache \
        /opt/conda/conda-bld \
        /opt/conda/pkgs \
        /home/mambauser/.cache \
        /home/mambauser/.cargo \
        /home/mambauser/.conda \
        /home/mambauser/.mamba \
        /home/mambauser/.rustup \
        /opt/capcruncher/.git \
        /opt/capcruncher/.pytest_cache \
        /opt/capcruncher/.ruff_cache \
        /opt/capcruncher/.uv-cache \
        /opt/capcruncher/build \
        /opt/capcruncher/tests \
        /opt/capcruncher/docs \
        /opt/capcruncher/dist \
        /tmp/* \
        /var/tmp/*

ENV CONDA_PREFIX=/opt/conda \
    CARGO_HOME=/tmp/cargo \
    MPLCONFIGDIR=/tmp/matplotlib \
    PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    PIP_CACHE_DIR=/tmp/pip-cache \
    PIP_NO_CACHE_DIR=1 \
    RUSTUP_HOME=/tmp/rustup \
    XDG_CACHE_HOME=/tmp/.cache \
    PYTHONUNBUFFERED=1

ENTRYPOINT ["capcruncher"]
CMD ["--help"]
