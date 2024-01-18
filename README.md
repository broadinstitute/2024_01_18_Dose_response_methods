# 2024_01_18_Dose_response_methods

Concentration-response analysis of omics data has emerged as a promising approach for integrating toxicogenomics data into a risk assessment context. In this approach, chemical exposures are conducted at multiple concentrations, omics data are measured from each sample, and curves are fit to each individual feature (i.e., a gene) to determine the concentration at which each feature surpasses levels observed in the control samples, also called the benchmark dose (BMD). The distribution of feature-level BMDs is then examined to identify the concentration at which the tested chemical causes a concerted perturbation across multiple features, called the omics point-of-departure (POD). While concentration-response analysis has been mainly applied to transcriptomics data, more recently it has been expanded to other omics types, including metabolomics and image-based cell profiling (Cell Painting) data. 

There are many different methods for deriving omics PODs, however, not all methods are suitable for all data modalities. In particular, popular methods that rely on molecular pathway annotations are especially difficult or impossible to generalize beyond transcriptomics data. Previous work has evaluated different methods for computing PODs with single omics data (e.g. transcriptional) by examining their reproducibility/stability and similarity to in vivo PODs derived from apical outcomes. The present work expands this investigation to a multi-omics context. There were two objectives: 1) to develop and evaluate a new method for POD calculation based on partial least squares (PLS), designed specifically to facilitate POD comparison across diverse omics data types, and 2) to evaluate the performance of both the novel and existing methods for POD calculation with multi-omics data.

## Documents
Include link to project folder once created. 

## What's in this repo?
This repo contains the analysis scripts for developing and evaluating dose-response analysis methods. The data is stored in a separate repo [PUT LINK ONCE CREATED], which is added as a submodule to this repo. The data are previously published transcriptomics and Cell Painting data, as described [here](https://www.sciencedirect.com/science/article/pii/S0041008X22001776). 

## How to use this repo?

1. Fork the repo
2. Clone the repo
  
    ``` bash
    git clone git@github.com:<YOUR USER NAME>/2024_01_18_Dose_response_methods.git
    ```
3. Download the contents of the submodule (Complete after initializing data repo)
