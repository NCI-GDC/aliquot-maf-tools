FROM python:3.5

# Copy over source
COPY . aliquot-maf-tools/
WORKDIR /aliquot-maf-tools

# Installing dependencies
RUN bash -c "./repo-install.sh" && \
    pip install -r docker-requirements.txt

# Install bio-submitter-qc
RUN pip install .
