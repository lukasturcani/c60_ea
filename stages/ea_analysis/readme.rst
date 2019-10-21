EA Analysis
===========

This stage is responsible for extracting the output of the EA into a
format suitable for analysis and then carrying out that analysis.
As such, the stage is divided into two sub-stages each found in the
``stages`` folder. The order in which the stages must be executed is

#. ``prepare_output`` - extracts output of the EA into suitable format
#. ``run_analysis`` - performs analysis, creates plots

Read the readme of each stage for further info.
