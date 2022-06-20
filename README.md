Supporting code for the article:

> AN Prybutok, JS Yu, JN Leonard\*, and N Bagheri\*. (2022). Mapping CAR T-cell design space using agent-based models. *Frontiers in Molecular Biosciences.* \*co-corresponding authorship.

## Summary of README contents

+ [Repository-overview](#repository-overview)
+ [Installation and running instructions](#installation-and-running-instructions)
  - [Install Python and requirements](#install-python-and-requirements)
  - [Clone GitHub repository](#clone-github-repository)
  - [Run Jupyter notebook](#run-jupyter-notebook)
+ [Pipeline summary](#pipeline-summary)
+ [`example/` directory contents](#`example/`-directory-contents)
  - [Directory overview](#directory-overview)
  - [`files/` directory contents](#`files/`-directory-contents)
    - [`toy/` directory contents](#`toy/`-directory-contents)
    - [`full/` directory contents](#`full/`-directory-contents)
  - [`figures/` directory contents](#`figures/`-directory-contents)


## Repository overview

The `scripts/` directory contains all the scripts used for analyzing the output of the simulations using [CARCADE v1.0](https://github.com/bagherilab/CARCADE).

The `pipeline.html` file walks through the data processing and plotting pipeline used to analyze the output of the simulations. Sample plots are included for each plotting section. This file shows the process in full and python commands used to generate all outputs. Go to https://bagherilab.github.io/carcade_mapping_design_space/pipeline.html to view the published html page.

The `pipeline_example.ipynb` file is a Jupyter Notebook that walks through an example of the data processing and plotting pipelines used to analyze the output of the simulations. In order to run this example notebook, clone this repository and install the packages listed in the `requirements.txt`, then run the various sections of interest in the notebook. No sections within the notebook are dependent on one another as we provide the outputs of all sections and the required folder structure used in the notebook.

The `examples/` directory contains data files and figures used in the `pipeline.html` and `pipeline_example.ipynb`. More detail about the contents of this directory are described below.

The `requirements.txt` file lists the package requirements for running this code. Additionally, code was run using Python version 3.6.8.

## Installation and running instructions

### Install Python and requirements

Install or create a virtual environment with Python version 3.6.8.

YOu can use `pip` to install the packages listed in the `requirements.txt` using the following command:

```
$ pip install -r requirements.txt
```

### Clone GitHub repository

In the command line, navigate to location on your computer where folder where you would like this repository stored and run the following command:

```
$ git clone https://github.com/bagherilab/carcade_mapping_design_space.git
```

### Run Jupyter notebook

To run the `pipeline_example.ipynb` jupyter notebook, navigate to the `carcade_mapping_design_space` repository folder and type the following command:

```
$ jupyter notebook
```

This will hopen up a web page with the contents of the repository. Click on the `pipeline_example` to run the notebook.

## Pipeline summary

The data processing pipeline progresses through the following stages, each of which has a corresponding folder with code required for this stage in the `scripts/` directory. Prior to entering this pipeline, output `.json` files are grouped by replicates (differing seeds of the same simulation setup) and compressed into `.tar.xz` files

+ **parse** - processes the CARCADE compressed output `.tar.xz` files into `.pkl` files with specific data structure for easy parsing
+ **analyze** - analyzes the parsed `.pkl` files for simulaton information over time for all cell types (outputs `.pkl` files); this stage is used to collect the following types of simulation information:
  - **cells** - cell counts, cell state counts/fractions, cell volumes, average cell cycle length, etc. This can be done for all cells (**cells**) or for only cells that share a location with at least one cancer cell (**sharedlocs**)
  - **environment** - molecule/nutrient concentrations over time
  - **spatial** - cell counts across simulation radius over time
  - **lysed** - tissue cell killing over time
+ **subset** - grabs subsets of analyzed `.pkl` simulation files within a given folder that match specified setup information and stores them in one combined `.pkl` file
+ **plot** - plots slices of data given subsetted `.pkl` file based on type of data contained
+ **stats** - analyzes whole data and outcomes of given subsetted `.pkl` file based on type of data contained
+ **image** - produces `.svg` images of cells or graph vasculature from given `.json` files

Each of the provided pipelines walks through these processes.

## `example/` directory contents

### Directory overview

The `examples/` directory contains the following folder structure:

  ```
  examples/
    |___files/
    |___figures/
  ````

where there are example data files in the `files/` folder and figure outputs from the pipeline in the `figures/` folder.

The `files/` folder contains example data files that are used in the pipeline. `toy` simulations are those that are short time scale (2-day) simulations that enable time efficient exploration of the processing part of the pipeline. `full` simulations are full-length simulations used in the paper to enable data plotting and analysis.

All figure outputs in the `figures/` folder are generated using the `full` data.

### `files/` directory contents

These example files are sorted into `full/` and `toy/` data files with the following structure:

```
files/
|___full/
|___toy/
```

The `toy` data files in the `toy/` directory are short time scale (2-day) simulations that enable showing how the data processing pipeline works in a reasonable time frame. These are used for the **parse**, **analyze**, and **subset** parts of the example pipeline.

The `full` data files in the `full/` directory contain model output `.json` files and subsetted data files directly used in the paper and as part of the **plot**, **stats**, and **image** part of the example pipeline.

Files within the `full/` and `toy/` folders are sorted into subfolders by context or data type, including `coculture/`, `tissue/`, and `ranked/`. Files in `coculture/` folders correspond to co-culture `dish` context simulations. Files in `tissue/` folder correspond to `tissue` context simulations. Files in the `ranked/` folder enable comparison between the performance rank of the effective treatments from the realistic co-culture `dish` when used in the `tissue` context.

#### `toy/` directory contents

`toy` simulations are those that are short time scale (2-day) simulations that enable time efficient exploration of the processing part of the pipeline.

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

The `tars/` folder contains the compressed `.json` file outputs from the CARCADE model, where files are zipped together with replicates of the same simulation setup.

The `parsed/` folder contains parsed `.pkl` files created from the `.tar` files.

The `lysis/`, folder contains `.json.LYSIS` file outputs from the CARCADE model.

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

#### `full/` directory contents

`full` simulations are full-length simulations used in the paper to enable data plotting and analysis.

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

The `subset/` folders contains the same structure and contents as the `subset/` folder in the toy file, except the contents are of full length simulations used in the paper.

The `jsons/` folders contain output `.json` files from the CARCADE model that are used in the imaging portion of the pipeline. The `coculture/` folder only contains `.json` files with cell location and type information. For the `tissue/` folder, these are separated into `cells/` and `graph/` `.json` files, where the `cells/` folder contains `.json` files with cell location and type information while the `graph/` folder contains `.json` files with vasculature graph information.

### `figures/` directory contents

The `figures/` folder contains `.svg` and `.pdf`  figure outputs from the analysis pipeline. These figures are sorted within this folder by type or context using the following folder structure:

```
figures/
|___coculture/
|___heuristics/
|___kill_curves/
|___ranked/
|___tissue/
```

The `coculture/` folder contains sorted figures corresponding to co-culture `dish` context simulations.

The `tissue/` folder contains sorted figures corresponding to correspond to `tissue` context simulations.

The `heuristics/` folder contains figure outputs describing the binding heuristics for CAR-antigen and PD1-PDL1 binding used in the model.

The `kill_curves/` folder contains figures showing the kill curves based on data from experimental literature.

The `ranked/` folder contains images comparing the rank of effective treatments from the realistic co-culture `dish` context to the `tissue` context.

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
