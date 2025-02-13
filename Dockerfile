FROM python:3.8-slim

# Metadata
LABEL software="DeepSig"
LABEL software.version="20250212"
LABEL description="an open source software tool to predict signal peptides in proteins"
LABEL website="https://deepsig.biocomp.unibo.it"
LABEL documentation="https://deepsig.biocomp.unibo.it"
LABEL license="GNU GENERAL PUBLIC LICENSE Version 3"
LABEL maintainer="Castrense Savojardo <castrense.savojardo2@unibo.it>"

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=true \
    TF_CPP_MIN_LOG_LEVEL=3 \
    DEEPSIG_ROOT=/usr/src/deepsig \
    PATH=/usr/src/deepsig:$PATH

WORKDIR /usr/src/deepsig
COPY . .

# Install DeepSig and its dependencies
RUN pip install --no-cache-dir .

# Use the installed DeepSig wrapper
ENTRYPOINT ["deepsig"]
