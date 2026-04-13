ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=4.4.1

FROM ${REGISTRY}/python3.12-builder:${BASE_CONTAINER_VERSION} as builder

# Create virtual environment for the builder
ENV VIRTUAL_ENV=/app

COPY ./ /aliquotmaf

WORKDIR /aliquotmaf

ARG PIP_INDEX_URL
RUN ls -la && git status && uv tool run --with tox-uv-bare tox -e build

# Install the built package into the virtual environment
RUN uv pip install --no-binary=pysam /aliquotmaf/dist/*.whl

FROM ${REGISTRY}/python3.12:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="aliquotmaf" \
      org.opencontainers.image.description="Tools for creating and filtering aliquot-level MAFs" \
      org.opencontainers.image.source="https://github.com/NCI-GDC/aliquot-maf-tools" \
      org.opencontainers.image.vendor="NCI GDC"

# Copy the virtual environment from builder (includes all installed packages)
COPY --from=builder /app /app

# Set up the virtual environment to be used
ENV VIRTUAL_ENV=/app \
    PATH="/app/bin:$PATH"

USER app
