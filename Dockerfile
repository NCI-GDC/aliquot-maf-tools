ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=4

FROM ${REGISTRY}/python3.12-builder:${BASE_CONTAINER_VERSION} AS builder

WORKDIR /app
RUN --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    dnf install htslib && \
    uv sync --no-install-project --no-dev --active

COPY . /app
RUN uv sync --no-dev --active

FROM ${REGISTRY}/python3.12-builder:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

RUN dnf install htslib
COPY --from=builder --chown=app /app /app

ENV PATH="/app/.venv/bin:$PATH"

USER app

RUN aliquotmaf --help

ENTRYPOINT ["/usr/bin/dumb-init", "--"]

CMD ["aliquotmaf", "--help"]
