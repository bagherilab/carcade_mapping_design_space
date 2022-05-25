Supporting code for the article:

> AN Prybutok, JS Yu, JN Leonard\*, and N Bagheri\*. (2022). Mapping CAR T-cell design space using agent-based models. *Frontiers in Molecular Biosciences.* \*co-corresponding authorship.

# Analysis scripts overview

- **[Repository overview](#repository-overview)**
- **[Pipeline walkthroughs](#pipeline-walkthroughs)**
  - [Full data processing and plotting pipeline](#full-data-processing-and-plotting-pipeline)
  - [Example data processing and plotting pipeline](#full-data-processing-and-plotting-pipeline)

## Repository overview

The `scripts/` directory contains all the scripts used for analyzing the output of the simulations using [CARCADE v1.0](#https://github.com/bagherilab/CARCADE)

The `pipeline_data_processing_and_plotting.html` and  `pipeline_data_processing_and_plotting_example.ipynb` files walk through the full data processing and plotting pipelines used to analyze the data from the paper.

The `examples/` directory contains example data files, sorted by type of data, used in the pipeline and figure outputs of the analysis scripts.

The `requirements.txt` file lists the package requirements for running this code. Additionally, code was run using Python version 3.6.8.

## Pipeline walkthroughs

### Full data processing and plotting pipeline

The `pipeline_data_processing_and_plotting.html` file walks through the data processing and plotting pipeline used to analyze the output of the simulations. Sample plots are included for each plotting section. This file is not runable, but is meant to show the process in full and python commands used to generate all outputs.

Go to https://bagherilab.github.io/carcade_mapping_design_space/pipeline_data_processing_and_plotting.html to view the published html page.

### Example data processing and plotting pipeline

The `pipeline_data_processing_and_plotting_example.ipynb` file is a runnable Jupyter Notebook that walks through an example of the data processing and plotting pipelines used to analyze the output of the simulations.
