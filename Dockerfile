FROM quay.io/ncigdc/bio-python:3.6

ENV BINARY=maflib

COPY ./dist/ /opt

WORKDIR /opt

RUN make init-pip \
  && ln -s /opt/bin/${BINARY} /bin/${BINARY} \
  && chmod +x /bin/${BINARY}

ENTRYPOINT ["/bin/aliquotmaf"]

CMD ["--help"]
