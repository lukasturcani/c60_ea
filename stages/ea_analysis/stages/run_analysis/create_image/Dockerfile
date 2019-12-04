FROM c60_ea_prepare_output:latest
WORKDIR /ea_analysis/run_analysis
RUN mkdir /ea_analysis/run_analysis/output
COPY . /ea_analysis/run_analysis
SHELL ["/bin/bash", "-c"]

RUN source activate c60_ea; \
    for N in {1..5}; \
    do \
    PYTHONPATH="/ea_analysis/run_analysis/code:/ea_analysis/code" \
    python /ea_analysis/run_analysis/code/create_building_block_plots.py \
    "/ea_analysis/ea_progress/$N.pkl" \
    "/ea_analysis/run_analysis/output/$N/plots"; \
    PYTHONPATH="/ea_analysis/run_analysis/code:/ea_analysis/code" \
    python /ea_analysis/run_analysis/code/generational_fittest.py \
    "/ea_analysis/ea_progress/$N.pkl" \
    "/ea_analysis/run_analysis/output/$N/best_structures"; \
    done; \
    PYTHONPATH="/ea_analysis/run_analysis/code:/ea_analysis/code" \
    python /ea_analysis/run_analysis/code/create_diversity_plot.py \
    "/ea_analysis/ea_progress/1.pkl" \
    "/ea_analysis/ea_progress/2.pkl" \
    "/ea_analysis/ea_progress/3.pkl" \
    "/ea_analysis/ea_progress/4.pkl" \
    "/ea_analysis/ea_progress/5.pkl" \
    "/ea_analysis/run_analysis/output/diversity_plot.png";
