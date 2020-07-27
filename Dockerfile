FROM continuumio/miniconda3:4.7.12
RUN pip install --upgrade pip && \
    pip install cirrocumulus==1.1.1
SHELL ["/bin/bash", "-c"]


