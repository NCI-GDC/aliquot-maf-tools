FROM docker.osdc.io/ncigdc/python3.8-builder as builder

COPY ./ /opt

WORKDIR /opt

RUN pip install tox && tox -p

FROM docker.osdc.io/ncigdc/python3.8

COPY --from=builder /opt/dist/*.tar.gz /opt
COPY requirements.txt /opt

WORKDIR /opt

RUN pip install -r requirements.txt \
	&& pip install *.tar.gz \
	&& rm -f *.tar.gz requirements.txt

ENTRYPOINT ["aliquot_maf_tools"]

CMD ["--help"]
