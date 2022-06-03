Supporting code for the article:

> AN Prybutok, JS Yu, JN Leonard\*, and N Bagheri\*. (2022). Mapping CAR T-cell design space using agent-based models. *Frontiers in Molecular Biosciences.* \*co-corresponding authorship.

## Repository overview

The `scripts/` directory contains all the scripts used for analyzing the output of the simulations using [CARCADE v1.0](https://github.com/bagherilab/CARCADE).

The `pipeline.html` file walks through the data processing and plotting pipeline used to analyze the output of the simulations. Sample plots are included for each plotting section. This file is not runable, but is meant to show the process in full and python commands used to generate all outputs. Go to https://bagherilab.github.io/carcade_mapping_design_space/pipeline.html to view the published html page.

The `pipeline_example.ipynb` file is a runnable Jupyter Notebook that walks through an example of the data processing and plotting pipelines used to analyze the output of the simulations. In order to run this example notebook, clone this repository and install the packages listed in the `requirements.txt`, then run the various sections of interest in the notebook. No sections within the notebook are dependent on one another as we provide the outputs of all sections and the required folder structure used in the notebook.

The `requirements.txt` file lists the package requirements for running this code. Additionally, code was run using Python version 3.6.8.

The `examples/` directory contains the following folder structure:

  ```
  examples/
    |___files/
    |     |___full/
    |     |___toy/
    |___figures/
  ````

The `examples/` directory contains example data files in the `files/` folder contains example data files that are used in the pipeline. These example files are sorted into `full/` data files, which are those directly used in the paper, and `toy/` data files, which are short time scale simulations that are not used in the paper but that are meant to show how the data processing pipeline works in a reasonable time frame.

Files within these folders are sorted into subfolders by context or data type, including `coculture/`, `tissue/`, and `ranked/`. Files in `coculture/` folders correspond to co-culture `dish` context simulations. Files in `tissue/` folder correspond to `tissue` context simulations. Files in the `ranked/` folder enable comparison between the performance rank of the effective treatments from the realistic co-culture `dish` when used in the `tissue` context.

The `toy/` directory contains the following folder structure:

```
toy/
|___coculture/
|     |___analyzed/
|     |___lysis/
|     |___parsed/
|     |___subset/
|     |___tars/
|___ranked/
|___tissue/
|     |___analyzed/
|     |___lysis/
|     |___parsed/
|     |___subset/
|     |___tars/
```

The `tars/` folder contains the compressed `.json` file outputs from the CARCADE model, where files are zipped together with replicates of the same simulation setup. The `parsed/` folder contains parsed `.pkl` files created from the `.tar` files. The `lysis/`, folder contains `.json.LYSIS` file outputs from the CARCADE model.

The `analyzed/` and `subset/` folders contain `.pkl` files, which are sorted into subfolders, corresponding to the cell analysis and subsetted analysis or lysis files that result from the analysis pipeline. The subfolder structure is as follows:

```
analyzed or subset/
|___cells/
|___environment/
|___lysed/
|___sharedlocs/
|___spatial/
````

Each of these folders contain outputs from the processing pipeline described in the pipeline walkthroughs.

The `full/` folder contains a similar structure to that of the `toy/` folder, as shown below:

```
full/
|___coculture/
|     |___jsons/
|     |___subset/
|___ranked/
|___tissue/
|     |___jsons/
|     |     |___cells/
|     |     |___graph/
|     |___subset/
```

The `subset/` folders contains the same structure and contents as the `subset/` folder in the toy file, except the contents are of full simulations used in the paper.

The `jsons/` folders contain output `.json` files from the CARCADE model that are used in the imaging portion of the pipeline. The `coculture/` folder only contains `.json` files with cell location and type information. For the `tissue/` folder, these are separated into `cells/` and `graph/` `.json` files, where the `cells/` folder contains `.json` files with cell location and type information while the `graph/` folder contains `.json` files with vasculature graph information.

The `figures/` folder contains `.svg` and `.pdf`  figure outputs from the analysis pipeline. These figures are sorted within this folder by type or context using the following folder structure:

```
figures/
|___coculture/
|___heuristics/
|___kill_curves/
|___ranked/
|___tissue/
```

The `coculture/` and `tissue/` folders contain the following subfolders:

```
coculture or tissue/
|___cells/
|___environment/
|___images/
|___lysed/
|___sharedlocs/
|___spatial/
|___stats/
```

Each of these folders contain figure outputs from the processing pipeline described in the pipeline walkthroughs.

The `heuristics/` folder contains figure outputs describing the binding heuristics for CAR-antigen and PD1-PDL1 binding used in the model. The `kill_curves/` folder contains figures showing the kill curves based on data from experimental literature. The `ranked/` folder contains images comparing the rank of effective treatments from the realistic co-culture `dish` context to the `tissue` context.
