ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=4

FROM ${REGISTRY}/python3.12-builder:${BASE_CONTAINER_VERSION} AS builder

ARG HTSLIB_VERSION=1.23.1
ADD https://nexus.osdc.io/repository/github/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 /htslib.tar.bz2
WORKDIR /
RUN tar xf htslib.tar.bz2 -C .
WORKDIR /htslib-${HTSLIB_VERSION}
RUN <<EOR
./configure --prefix=/usr/local
make
make install
EOR

ARG SAMTOOLS_VERSION=1.23.1
ADD https://nexus.osdc.io/repository/github/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 /samtools.tar.bz2
WORKDIR /
RUN tar xf samtools.tar.bz2 -C .
WORKDIR /samtools-${SAMTOOLS_VERSION}
RUN <<EOR
./configure --prefix=/usr/local
make
make install
EOR

ENV HTSLIB_LIBRARY_DIR=/usr/local/lib \
    HTSLIB_INCLUDE_DIR=/usr/local/include

WORKDIR /app
RUN --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --no-install-project --no-dev --active

COPY . /app
RUN uv sync --no-dev --active

FROM ${REGISTRY}/python3.12:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

COPY --from=builder /usr/local /usr/local
COPY --from=builder --chown=app /app /app

ENV PATH="/app/.venv/bin:$PATH"

USER app

RUN aliquotmaf --help

ENTRYPOINT ["/usr/bin/dumb-init", "--"]

CMD ["aliquotmaf", "--help"]
