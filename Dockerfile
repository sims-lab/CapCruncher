FROM mambaorg/micromamba:2.0.5

ARG CAPCRUNCHER_VERSION=0.0.0+container
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG QUARTO_VERSION=1.9.37

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 \
        pip \
        'bedtools>=2.31.0' \
        'bowtie2>=2.4.4' \
        coreutils \
        curl \
        cxx-compiler \
        'fastqc<=0.12.1' \
        flash2 \
        git \
        'iced>=0.5.10' \
        jupyterlab \
        pairix \
        pigz \
        rust \
        'samtools>=1.15.1' \
        'trim-galore<=0.6.10' \
        ucsc-bedgraphtobigwig \
        ucsc-bedtobigbed && \
    ln -sf /opt/conda/bin/flash2 /opt/conda/bin/flash && \
    micromamba clean --all --yes

WORKDIR /opt/capcruncher

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/capcruncher

RUN printf 'setuptools<80\n' > /tmp/pip-build-constraints.txt && \
    python -m pip install --no-cache-dir --upgrade pip && \
    PIP_CONSTRAINT=/tmp/pip-build-constraints.txt \
    SETUPTOOLS_SCM_PRETEND_VERSION_FOR_CAPCRUNCHER="${CAPCRUNCHER_VERSION}" \
        python -m pip install --no-cache-dir '.[full]'

RUN QUARTO_ARCH="$(uname -m)" && \
    case "${QUARTO_ARCH}" in \
        aarch64|arm64) QUARTO_ARCH="arm64" ;; \
        x86_64|amd64) QUARTO_ARCH="amd64" ;; \
        *) echo "Unsupported Quarto architecture: ${QUARTO_ARCH}" >&2; exit 1 ;; \
    esac && \
    mkdir -p "/opt/conda/share/quarto/${QUARTO_VERSION}" && \
    curl -fsSL \
        "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${QUARTO_ARCH}.tar.gz" \
        -o /tmp/quarto.tar.gz && \
    tar -xzf /tmp/quarto.tar.gz \
        -C "/opt/conda/share/quarto/${QUARTO_VERSION}" \
        --strip-components=1 && \
    ln -sf "/opt/conda/share/quarto/${QUARTO_VERSION}/bin/quarto" /opt/conda/bin/quarto && \
    rm /tmp/quarto.tar.gz && \
    quarto --version

ENV CONDA_PREFIX=/opt/conda \
    MPLCONFIGDIR=/tmp/matplotlib \
    PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    XDG_CACHE_HOME=/tmp/.cache \
    PYTHONUNBUFFERED=1

ENTRYPOINT ["capcruncher"]
CMD ["--help"]
