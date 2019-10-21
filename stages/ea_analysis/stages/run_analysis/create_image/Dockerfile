FROM c60_ea_prepare_output:latest
WORKDIR /ea_analysis/run_analysis
RUN mkdir /ea_analysis/run_analysis/output
COPY . /ea_analysis/run_analysis
SHELL ["/bin/bash", "-c"]

RUN source activate c60_ea; \
    for N in {1..5}; \
    do \
    PYTHONPATH="/ea_analysis/run_analysis/code" \
    python /ea_analysis/run_analysis/code/create_building_block_plots.py \
    "/ea_analysis/ea_progress/$N.pkl" \
    "/ea_analysis/run_analysis/output/$N"; \
    done
