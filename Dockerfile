# Multi-stage Dockerfile for qchem-vqe (Linux / macOS Docker Desktop / HPC Apptainer)
FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake git ca-certificates gfortran \
    libopenblas-dev liblapacke-dev libnlopt-dev libcjson-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .
RUN make fetch-libcint && make all && make test && make install PREFIX=/opt/qchem-vqe

FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates libgomp1 libgfortran5 libquadmath0 \
    libopenblas0 liblapacke libnlopt0 libcjson1 libcurl4t64 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/qchem-vqe /usr/local
ENV LD_LIBRARY_PATH="/usr/local/lib"
ENV QCHEM_DATA_DIR="/usr/local/share/qchem-vqe/data"
WORKDIR /workspace
ENTRYPOINT ["qchem-vqe"]
