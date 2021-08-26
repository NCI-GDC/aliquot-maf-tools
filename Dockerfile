FROM quay.io/ncigdc/python36-builder as builder

COPY ./ /opt

WORKDIR /opt

RUN pip install tox && tox -p

FROM quay.io/ncigdc/python36

COPY --from=builder /opt/dist/*.tar.gz /opt
COPY requirements.txt /opt

WORKDIR /opt

RUN pip install -r requirements.txt \
	&& pip install *.tar.gz \
	&& rm -f *.tar.gz requirements.txt

ENTRYPOINT ["aliquotmaf"]

CMD ["--help"]
