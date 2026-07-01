# DNA-origami-SMFS-CaE-analysis

This repository contains custom analysis code associated with the protocol:

**DNA origami-enhanced single-molecule force spectroscopy and AlphaFold analysis of protein folding**

The code supports optical-tweezers trace processing, passive-mode segmentation, hidden Markov model (HMM) fitting, dwell-time analysis, contour-length calculation, AF-Cluster conformational-ensemble analysis, and AF-Cluster post-processing for Cα-based RMSD, pairwise-distance PCA and figure-source-data generation.

## Repository contents

```text
DNA-origami-SMFS-CaE-analysis/
├── README.md
├── LICENSE
├── environment.yml
├── requirements.txt
├── matlab/
│   └── passive_hmm/
│       ├── LA_hires_passive_processing.m
│       ├── LA_hires_passive_split.m
│       ├── isomerization_split.m
│       ├── LA_pyHMM_oneMolecule.m
│       └── supporting_functions/
├── notebooks/
│   └── AFcluster_CC.ipynb
├── python/
│   └── afcluster_postprocess.py
├── docs/
│   └── software_usage_instructions.md
├── example_data/
│   └── README.md
└── legacy/
    └── LA_pyHMM_oneConstruct.m
```

## Code modules

### 1. Passive-mode HMM analysis

The MATLAB scripts in `matlab/passive_hmm/` are used for selecting passive-mode regions, splitting kinetically distinct portions of single-molecule trajectories when necessary, fitting molecule-level HMMs through `hmmlearn`, and exporting state assignments and dwell times.

Main scripts:

- `LA_hires_passive_processing.m`
- `LA_hires_passive_split.m`
- `isomerization_split.m`
- `LA_pyHMM_oneMolecule.m`

The script `LA_pyHMM_oneConstruct.m` is provided only as a legacy script and is not used for the DNA origami-handle datasets because its implementation assumes WLC handles.

### 2. AF-Cluster conformational-ensemble analysis

The notebook `notebooks/AFcluster_CC.ipynb` is a modified AF-Cluster notebook used for DBSCAN-based MSA clustering and AlphaFold2 ensemble generation.

The analysis used in the protocol applies:

- the internal >25% MSA gap filter implemented in the AF-Cluster workflow;
- AlphaFold2 `model_3_ptm`;
- three recycles;
- random seed 0.

### 3. AF-Cluster post-processing

The Python script `python/afcluster_postprocess.py` is used to calculate:

- mean pLDDT values;
- Cα RMSDs relative to a common reference structure;
- all-pair Cα–Cα distance features;
- principal-component coordinates;
- machine-readable output tables;
- figure-source-data files and plots.

## Requirements

The analysis was developed and tested using the following software environment:

- MATLAB R2018a or later
- Python 3.8 or later
- hmmlearn 0.2.7 or later
- ColabFold 1.5.2 or later
- MMseqs2
- MDTraj 1.9.9 or later
- scikit-learn 1.2.2 or later
- ChimeraX 1.7 or later

Python dependencies are listed in `requirements.txt`. A conda environment file is provided as `environment.yml`.

## Installation

Clone the repository:

    git clone https://github.com/JingyiYe/DNA-origami-SMFS-CaE-analysis.git
    cd DNA-origami-SMFS-CaE-analysis

Create the Python environment:

    conda env create -f environment.yml
    conda activate smfs-cae-analysis

Alternatively, install the Python dependencies using:

    pip install -r requirements.txt

MATLAB scripts can be run directly after adding `matlab/passive_hmm/` and its supporting functions to the MATLAB path.

## Input files

Expected input files include:

1. passive-mode optical-tweezers force-time traces exported from the acquisition software;
2. segmented passive-mode regions for HMM analysis;
3. AF-Cluster-generated PDB files;
4. AlphaFold2 or AF-Cluster prediction metadata;
5. a reference PDB structure for Cα RMSD calculation.

Detailed input-file formats are described in `docs/software_usage_instructions.md`.

## Output files

The scripts generate:

1. HMM state assignments;
2. dwell-time distributions;
3. transition-rate tables;
4. contour-length analysis tables;
5. Cα RMSD tables;
6. pairwise Cα–Cα distance matrices;
7. PCA coordinates;
8. figure-source-data tables;
9. diagnostic plots.

## Runnable example

A representative runnable example, including example input files, expected outputs and execution instructions, is deposited in figshare:

**figshare DOI:** TODO_FIGSHARE_DOI

The example is intended to demonstrate the analysis workflow and file formats. It does not include the complete raw optical-tweezers dataset.

## Version

The analyses described in the protocol correspond to the following GitHub commit:

**Commit:** TODO_FULL_COMMIT_HASH

## License

This repository is released under the **MIT License**.

Third-party software, including AF-Cluster, ColabFold, MMseqs2, hmmlearn, MDTraj, scikit-learn, ChimeraX and MATLAB, is available from the respective developers and is subject to its original license.

## Associated manuscript

This repository is associated with an unpublished protocol manuscript:

**DNA origami-enhanced single-molecule force spectroscopy and AlphaFold analysis of protein folding**

A full citation will be added after publication.

## Contact

For questions about the code or the associated protocol, please contact:

**Honglu Zhang**  
Shanghai Jiao Tong University  
[z.hl@sjtu.edu.cn](mailto:z.hl@sjtu.edu.cn)
