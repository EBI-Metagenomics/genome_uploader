FROM mambaorg/micromamba:2.5.0

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER *.whl /tmp/

RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba clean -a -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN pip install --no-cache-dir /tmp/*.whl

ENTRYPOINT ["genome_upload"]