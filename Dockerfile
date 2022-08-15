ARG REGISTRY=quay.io
ARG BASE_CONTAINER_VERSION=2.0.1

FROM ${REGISTRY}/ncigdc/python3.10-builder:${BASE_CONTAINER_VERSION} as builder

COPY ./ /opt/build

WORKDIR /opt/build

RUN pip install tox && tox -e build

FROM ${REGISTRY}/ncigdc/python3.10-builder:${BASE_CONTAINER_VERSION}

COPY --from=builder /opt/build/dist/*.tar.gz /opt/build/
COPY requirements.txt /opt/build/

WORKDIR /opt/build

ARG PIP_EXTRA_INDEX_URL

RUN --mount=type=ssh pip install --no-deps -r requirements.txt \
	&& pip install --no-deps *.tar.gz \
	&& rm -f *.tar.gz requirements.txt

ENTRYPOINT ["aliquot_maf_tools"]

CMD ["--help"]
