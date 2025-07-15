ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM ${REGISTRY}/python3.12-builder:feat_uv-base-container as builder

WORKDIR /app
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-dev --active

COPY . /app
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --active

FROM ${REGISTRY}/python3.12:feat_uv-base-container

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

COPY --from=builder /venv /venv

USER app
