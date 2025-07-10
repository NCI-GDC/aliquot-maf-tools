ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM ${REGISTRY}/python3.12-builder:${BASE_CONTAINER_VERSION} as builder

COPY --from=ghcr.io/astral-sh/uv:0.7.12 /uv /uvx /bin/

COPY ./ /aliquotmaf

WORKDIR /aliquotmaf

RUN uv build --no-binary

WORKDIR /deps

RUN uv venv -p python3.12 && source .venv/bin/activate \
      && uvx pip wheel --no-binary pysam -r /aliquotmaf/requirements.txt

FROM ${REGISTRY}/python3.12:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

COPY --from=ghcr.io/astral-sh/uv:0.7.12 /uv /uvx /bin/
COPY --from=builder /aliquotmaf/dist/*.whl /aliquotmaf/
COPY --from=builder /deps/*.whl /aliquotmaf

WORKDIR /aliquotmaf

RUN uv pip install --no-deps *.whl \
	&& rm -f *.whl

USER app
