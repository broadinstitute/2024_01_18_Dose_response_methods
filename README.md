# 2024_01_18_Dose_response_methods

Concentration-response analysis of omics data has emerged as a promising approach for integrating toxicogenomics data into a risk assessment context. In this approach, chemical exposures are conducted at multiple concentrations, omics data are measured from each sample, and curves are fit to each individual feature (i.e., a gene) to determine the concentration at which each feature surpasses levels observed in the control samples, also called the benchmark dose (BMD). The distribution of feature-level BMDs is then examined to identify the concentration at which the tested chemical causes a concerted perturbation across multiple features, called the omics point-of-departure (POD). While concentration-response analysis has been mainly applied to transcriptomics data, more recently it has been expanded to other omics types, including metabolomics and image-based cell profiling (Cell Painting) data. 

There are many different methods for deriving omics PODs, however, not all methods are suitable for all data modalities. In particular, popular methods that rely on molecular pathway annotations are especially difficult or impossible to generalize beyond transcriptomics data. Previous work has evaluated different methods for computing PODs with single omics data (e.g. transcriptional) by examining their reproducibility/stability and similarity to in vivo PODs derived from apical outcomes. The present work expands this investigation to a multi-omics context. There were two objectives: 1) to develop and evaluate a new method for POD calculation based on partial least squares (PLS), designed specifically to facilitate POD comparison across diverse omics data types, and 2) to evaluate the performance of both the novel and existing methods for POD calculation with multi-omics data.

## Documents
All project documents are in the [this Google Drive folder](https://drive.google.com/drive/folders/1HY9FMRH9JYs3KQoDViD-ZLrlZhM25z6A?usp=drive_link).  

## What's in this repo?
This repo contains the analysis scripts for developing and evaluating dose-response analysis methods. The data is stored in a separate repo [`2024_01_18_Dose_response_methods-data`](https://github.com/broadinstitute/2024_01_18_Dose_response_methods-data), which is added as a submodule to this repo. The data are previously published transcriptomics and Cell Painting data, as described [here](https://www.sciencedirect.com/science/article/pii/S0041008X22001776). 

## How to use this repo for the first time?

1. Fork the repo
2. Clone the repo
  
    ``` bash
    git clone git@github.com:<YOUR USER NAME>/2024_01_18_Dose_response_methods.git
    ```
3. Download the contents of the submodule. If these steps give an "Unable to locate credentials" error, it probably means that you need to install the AWS CLI and configure with your credentials.

    ```bash
    git submodule update --init --recursive
    cd 2024_01_18_Dose_response_methods-data
    dvc pull
    ```  

## How to update this repo?
This repo is structured as a parent (2024_01_18_Dose_response_methods) with a submodule containing the data (2024_01_18_Dose_response_methods-data), which is tracked with dvc. All three components (parent repo, submodule repo, actual data files) must by synchronized. Here's how:

1. Navigate to 2024_01_18_Dose_response_methods and pull most recent changes:
    ``` bash
    git pull
    ```
2. Navigate to the submodule and pull most recent changes in both the repository and actual data files:
    ``` bash
    cd ./2024_01_18_Dose_response_methods-data
    git pull
    dvc pull
    ```

## How to sync modifications to the data files?
If you've changed the data files, updates must be committed and pushed to both S3 via dvc, to the relevant .dvc file within the data submodule repository, and to the parent repository such that it is linked to the most recent submodule commit. Here's how to fully update all three components after modifying the data files:

1. Navigate to the data submodule within the parent repository
2. Stage and push changes to the data to S3 using dvc. Since our main directory containing all dvc tracked data files is called 'data', the commands are:
    ``` bash
    dvc add data
    dvc push
    ```
3. Stage, commit, and push changes to the data.dvc file within the data submodule repository:
    ``` bash
    git add data.dvc
    git commit -am "commit message"
    git push
    ```
4. Navigate to the parent repository and stage, commit, and push changes to the submodule:
    ``` bash
    git add 2024_01_18_Dose_response_methods-data
    git commit -am "update submodule to most recent commit"
    git push
    ```
