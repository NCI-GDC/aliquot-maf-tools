ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM ${REGISTRY}/python3.12-builder:${BASE_CONTAINER_VERSION} as builder

COPY ./ /aliquotmaf

WORKDIR /aliquotmaf

ARG PIP_INDEX_URL
RUN pip install tox && tox -e build

WORKDIR /deps

ENV PIP_NO_BINARY=pysam

RUN pip wheel -r /aliquotmaf/requirements.txt

FROM ${REGISTRY}/python3.12:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

COPY --from=builder /aliquotmaf/dist/*.whl /aliquotmaf/
COPY --from=builder /deps/*.whl /aliquotmaf

WORKDIR /aliquotmaf

RUN pip install --no-deps *.whl \
	&& rm -f *.whl

USER app
