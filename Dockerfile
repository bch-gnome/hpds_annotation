FROM ubuntu:22.04

RUN apt-get update && apt-get install -y bcftools python3-minimal python3-pip

RUN pip install pysam

COPY transform_csq.v3.py /
