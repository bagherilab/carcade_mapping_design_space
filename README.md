Supporting code for the article:

> AN Prybutok, JS Yu, JN Leonard\*, and N Bagheri\*. (2022). Mapping CAR T-cell design space using agent-based models. *Frontiers in Molecular Biosciences.* 9:849363. doi: 10.3389/fmolb.2022.849363 \*co-corresponding authorship.

## Summary of README contents

+ [Repository overview](#repository-overview)
+ [Installation and running instructions](#installation-and-running-instructions)
  - [Install Python and requirements](#install-python-and-requirements)
  - [Clone GitHub repository](#clone-github-repository)
  - [Run Jupyter notebook](#run-jupyter-notebook)
+ [Pipeline summary](#pipeline-summary)
+ [`examples/` directory contents](#examples-directory-contents)
  - [Directory overview](#directory-overview)
  - [`files/` directory contents](#files-directory-contents)
    - [`toy/` directory contents](#toy-directory-contents)
    - [`full/` directory contents](#full-directory-contents)
  - [`figures/` directory contents](#figures-directory-contents)


## Repository overview

The `scripts/` directory contains all the scripts used for analyzing the output of the simulations using [CARCADE v1.0](https://github.com/bagherilab/CARCADE).

The `pipeline.html` file walks through the data processing and plotting pipeline used to analyze the output of the simulations. Sample plots are included for each plotting section. This file shows the process in full and Python functions and commands used to generate all outputs. Go to https://bagherilab.github.io/carcade_mapping_design_space/pipeline.html to view the published html page.

The `pipeline_example.ipynb` file is a Jupyter notebook that walks through examples of the main components of the data processing and plotting pipelines used to analyze the output of the simulations.

The `examples/` directory contains data files and figures used in the `pipeline.html` and `pipeline_example.ipynb`.

The `requirements.txt` file lists the package requirements for running this code.

## Installation and running instructions

### Install Python and requirements

Install or create a virtual environment with Python version 3.6.8.

You can use `pip` to install the packages listed in the `requirements.txt` using the following command:

```
$ pip install -r requirements.txt
```

### Clone GitHub repository

In the command line, navigate to the location on your computer where you would like this repository stored and run the following command:

```
$ git clone https://github.com/bagherilab/carcade_mapping_design_space.git
```

### Run Jupyter notebook

To run the `pipeline_example.ipynb` Jupyter notebook, navigate to the `carcade_mapping_design_space` repository folder and type the following command:

```
$ jupyter notebook
```

This will open up a web page with the contents of the repository where you can click on the `pipeline_example.ipynb` file to open the notebook.

No sections within the notebook are dependent on one another as we provide the outputs of all sections and the required folder structure used in the notebook.

## Pipeline summary

The data processing pipeline progresses through the following stages, each of which has a corresponding folder with code required for this stage in the `scripts/` directory. Prior to entering this pipeline, CARCADE model output simulation `.json` files from the growth profiler are grouped by replicates (differing seeds of the same simulation setup) and compressed into `.tar.xz` files

+ **parse** - processes the CARCADE compressed output `.tar.xz` files into `.pkl` files with a specific data structure for easy parsing
+ **analyze** - analyzes the parsed `.pkl` files for simulation information over time for all cell types (outputs `.pkl` files); this stage is used to collect the following types of simulation information:
  - **cells** - cell counts, cell state counts/fractions, cell volumes, average cell cycle length, etc. This can be done for all cells (**cells**) or for only cells that share a location with at least one cancer cell (**sharedlocs**)
  - **environment** - molecule/nutrient concentrations over time
  - **spatial** - cell counts across simulation radius over time
  - **lysed** - tissue cell killing over time
+ **subset** - grabs subsets of analyzed `.pkl` simulation files within a given folder that match specified setup information and stores them in one combined `.pkl` file
+ **plot** - plots slices of data given subsetted `.pkl` file based on type of data contained
+ **stats** - analyzes and plots whole data and outcomes of given subsetted `.pkl` file based on type of data contained
+ **image** - produces `.svg` images of tissue cells or graph vasculature at specified time points from given `.json` files

Each of the provided pipelines walks through these processes.

## `examples/` directory contents

### Directory overview

The `examples/` directory contains the following folder structure:

  ```
  examples/
    |___files/
    |___figures/
  ````

where there are example data files used as inputs for the pipeline are stored in the `files/` folder and figure outputs from the pipeline are stored in the `figures/` folder.

+ The `files/` folder contains two types of example data files that are used in the pipeline.
  - `toy` simulations are those that are short time scale (2-day) simulations that enable time-efficient exploration of the processing part of the pipeline. These are used for the **parse**, **analyze**, and **subset** parts of the example pipeline.
  - `full` simulations are full-length simulations used in the paper to enable data plotting and analysis. These include model output `.json` files and subsetted data files directly used in the paper and as part of the **plot**, **stats**, and **image** part of the example pipeline.
+ All figure outputs in the `figures/` folder are generated using the `full` data.

### `files/` directory contents

`full` and `toy` example files are sorted into `full/` and `toy/` directories, respectively, within the `files/` directory with the following structure:

```
files/
|___full/
|___toy/
```

where files within the `full/` and `toy/` directories are sorted into subfolders by context or data type, including `coculture/`, `tissue/`, and `ranked/`.

+ Files in `coculture/` folders correspond to co-culture `dish` context simulations.
+ Files in the `tissue/` folder correspond to `tissue` context simulations.
+ Files in the `ranked/` folder enable comparison between the performance rank of the effective treatments from the realistic co-culture `dish` when used in the `tissue` context.

#### `toy/` directory contents

The `toy/` directory contains the following folder structure:

```
toy/
|___coculture/
|     |___analyzed/
|     |___lysis/
|     |___parsed/
|     |___subset/
|     |___tars/
```

where each of the subfolders contains `toy` co-culture `dish` context simulation files for different aspects of the example pipeline:

+ The `tars/` folder contains the compressed `.json` file outputs from the CARCADE model growth profiler.
+ The `parsed/` folder contains parsed `.pkl` files created from the `.tar` files.
+ The `lysis/`, folder contains `.LYSIS.json` file outputs from the CARCADE model lysis profiler.
+ The `analyzed/` and `subset/` folders contain `.pkl` files, which are sorted into subfolders, corresponding to the cell analysis and subsetted analysis or lysis files that result from the analysis pipeline.

The subfolder structure for the `analyzed/` or `subset/` folders is as follows:

```
analyzed or subset/
|___cells/
|___environment/
|___lysed/
|___sharedlocs/
|___spatial/
````

#### `full/` directory contents

The `full/` folder contains the following folder structure:

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

where the subfolders within the `coculture/` and `tissue/` contains `full` co-culture `dish` or `tissue` context files, respectively, for different aspects of the example pipeline:

+ The `subset/` folders contains the same structure and contents as the `toy/coculture/subset/` folder, except the contents are of full-length simulations used in the paper.
+ The `jsons/` folders contain CARCADE output `.json` files that are used in the imaging portion of the pipeline.
  - The `coculture/` folder only contains `.json` files of co-culture `dish` context simulations with cell location and type information.
  - For the `tissue/` folder, these are separated into `cells/` and `graph/` `.json` files for `tissue` context simulations, where the `cells/` folder contains `.json` files with cell location and type information and result from the growth profiler while the `graph/` folder contains `.GRAPH.json` files with vasculature graph information and result from the graph profiler.

Figures generated from these `full` files are shown in the `pipeline.html`.

### `figures/` directory contents

The `figures/` folder contains `.svg` and `.pdf` figure outputs from the analysis pipeline. These figures are sorted within this folder by type or context using the following folder structure:

```
figures/
|___coculture/
|___heuristics/
|___kill_curves/
|___ranked/
|___tissue/
```

where these type or context subfolders contain the following:

+ The `coculture/` folder contains sorted figures corresponding to co-culture `dish` context simulations.
+ The `tissue/` folder contains sorted figures corresponding to correspond to `tissue` context simulations.
+ The `heuristics/` folder contains figure outputs describing the binding heuristics for CAR-antigen and PD1-PDL1 binding used in the model.
+ The `kill_curves/` folder contains figures showing the kill curves based on data from experimental literature.
+ The `ranked/` folder contains images comparing the rank of effective treatments from the realistic co-culture `dish` context to the `tissue` context.

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
