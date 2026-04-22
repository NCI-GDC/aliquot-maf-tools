ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=4

FROM ${REGISTRY}/amzn2023-builder:${BASE_CONTAINER_VERSION}

ENV UV_PYTHON=3.12

WORKDIR /app
RUN --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --no-install-project --no-dev --active --no-binary

COPY . /app
RUN uv sync --no-dev --active --no-binary

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

RUN chown -R app:app /app

ENV PATH="/app/.venv/bin:$PATH"

USER app

RUN aliquotmaf --help

ENTRYPOINT ["/usr/bin/dumb-init", "--"]

CMD ["aliquotmaf", "--help"]
