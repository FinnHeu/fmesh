FROM mambaorg/micromamba
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp
RUN sed -i 's/name: fmesh/name: base/g' /tmp/requirements.yml
RUN micromamba install -f /tmp/requirements.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1
CMD ["/tmp/fmesh.py"]
ENTRYPOINT ["/opt/conda/bin/python"]


