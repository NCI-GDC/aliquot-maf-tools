ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM ${REGISTRY}/python3.8-builder:${BASE_CONTAINER_VERSION} as builder

COPY ./ /aliquotmaf

WORKDIR /aliquotmaf

RUN pip install tox && tox -e build

FROM ${REGISTRY}/python3.8:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

COPY --from=builder /aliquotmaf/dist/*.whl /aliquotmaf/
COPY requirements.txt /aliquotmaf/

WORKDIR /aliquotmaf

RUN pip install --no-deps -r requirements.txt \
	&& pip install --no-deps *.whl \
	&& rm -f *.whl requirements.txt

USER app
