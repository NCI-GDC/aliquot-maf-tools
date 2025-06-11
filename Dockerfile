ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM ${REGISTRY}/python3.12-builder:${BASE_CONTAINER_VERSION} as builder

COPY ./ /aliquotmaf

WORKDIR /aliquotmaf

ENV UV_SYSTEM_PYTHON=1 UV_COMPILE_BYTECODE=1 UV_NO_MANAGED_PYTHON=1

RUN --mount=from=ghcr.io/astral-sh/uv,source=/uv,target=/bin/uv \
  uv tool run setuptools_scm \
  && uv build

FROM ${REGISTRY}/python3.12:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

WORKDIR /aliquotmaf

ENV UV_SYSTEM_PYTHON=1 UV_COMPILE_BYTECODE=1 UV_NO_MANAGED_PYTHON=1

RUN --mount=from=builder,source=/aliquotmaf/dist/,target=/aliquotmaf/dist/ \
    --mount=source=uv.lock,target=/aliquotmaf/uv.lock \
    --mount=source=pyproject.toml,target=/aliquotmaf/pyproject.toml \
    --mount=from=ghcr.io/astral-sh/uv,source=/uv,target=/bin/uv \
      uv sync --locked --no-dev --active \
	&& uv pip install --no-deps ./dist/*.whl

RUN ls -l \
  && ls -l ./*

USER app
