ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM ${REGISTRY}/python3.12-builder:${BASE_CONTAINER_VERSION} as builder

COPY ./ /aliquotmaf

WORKDIR /aliquotmaf

ENV UV_SYSTEM_PYTHON=1 UV_COMPILE_BYTECODE=1 UV_NO_MANAGED_PYTHON=1

RUN --mount=from=ghcr.io/astral-sh/uv,source=/uv,target=/bin/uv \
  pip install setuptools_scm \
  && python -m setuptools_scm \
  && uv build

FROM ${REGISTRY}/python3.12:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

WORKDIR /aliquotmaf

RUN --mount=from=builder,source=/aliquotmaf/dist/,target=dist/ \
    --mount=source=requirements.txt,target=requirements.txt \
      pip install --no-deps -r requirements.txt \
      pip install --no-deps ./dist/*.whl

USER app
