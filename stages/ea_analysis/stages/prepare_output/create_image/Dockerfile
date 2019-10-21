FROM continuumio/miniconda3:4.7.10
WORKDIR /ea_analysis
COPY . /ea_analysis
SHELL ["/bin/bash", "-c"]
RUN mkdir /ea_analysis/ea_progress
RUN conda env create -f environment.yaml

RUN source activate c60_ea; \
    for N in {1..5}; \
    do \
    PYTHONPATH="/ea_analysis/code" \
    python /ea_analysis/code/create_progress.py \
    "/ea_analysis/ea_output/$N/" \
    "/ea_analysis/ea_progress/$N.pkl"; \
    done
