FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml python=3.12 pip && \
    micromamba clean --all --yes

WORKDIR /opt/capcruncher

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/capcruncher

RUN python -m pip install --upgrade pip && \
    python -m pip install '.[full]'

ENV MPLCONFIGDIR=/tmp/matplotlib \
    XDG_CACHE_HOME=/tmp/.cache \
    PYTHONUNBUFFERED=1

ENTRYPOINT ["capcruncher"]
