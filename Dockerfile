FROM ubuntu:18.04

RUN apt-get update && apt-get install -y bcftools python-minimal

COPY transform_csq.v2.py /
