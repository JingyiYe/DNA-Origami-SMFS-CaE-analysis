# DNA-Origami-SMFS-CaE-analysis

This repository contains custom analysis code associated with the protocol:

**DNA origami-enhanced single-molecule force spectroscopy and AlphaFold analysis of protein folding**

The code supports optical-tweezers trace processing, passive-mode segmentation, hidden Markov model (HMM) fitting, dwell-time analysis, contour-length calculation, AF-Cluster conformational-ensemble analysis, and AF-Cluster post-processing for Cα-based RMSD, pairwise-distance PCA and figure-source-data generation.

## Repository contents

```text
DNA-Origami-SMFS-CaE-analysis/
├── README.md
├── LICENSE
├── environment.yml
├── requirements.txt
├── matlab/
│   └── passive_hmm/
│       ├── README.md
│       ├── hires_passive_processing.m
│       ├── hires_passive_split.m
│       ├── isomerization_split.m
│       └── pyHMM_oneMolecule.m
├── AF-Cluster_notebook/
│   ├── README.md
│   └── AFcluster.ipynb
└── python/
    ├── README.md
    └── AFcluster_postprocess.py
```

## Code modules

### 1. Passive-mode HMM analysis

The MATLAB scripts in `matlab/passive_hmm/` are used for selecting passive-mode regions, splitting kinetically distinct portions of single-molecule trajectories when necessary, fitting molecule-level HMMs through `hmmlearn`, and exporting state assignments and dwell times.

Main scripts:

- `hires_passive_processing.m`
- `hires_passive_split.m`
- `isomerization_split.m`
- `pyHMM_oneMolecule.m`

The script `pyHMM_oneConstruct.m` is provided only as a legacy script and is not used for the DNA origami-handle datasets because its implementation assumes WLC handles.

### 2. AF-Cluster conformational-ensemble analysis

The notebook `AF-Cluster_notebook/AFcluster.ipynb` is a modified AF-Cluster notebook used for DBSCAN-based MSA clustering and AlphaFold2 ensemble generation.

The analysis used in the protocol applies:

- the internal >25% MSA gap filter implemented in the AF-Cluster workflow;
- AlphaFold2 `model_3_ptm`;
- three recycles;
- random seed 0.

### 3. AF-Cluster post-processing

The Python script `python/AFcluster_postprocess.py` is used to calculate:

- mean pLDDT values;
- Cα RMSDs relative to a common reference structure;
- all-pair Cα–Cα distance features;
- principal-component coordinates;
- machine-readable output tables;
- figure-source-data files and plots.

## Requirements

The analysis was developed and tested using the following software environment:

- MATLAB R2018a
- Python 3.8
- hmmlearn 0.2.7
- ColabFold 1.5.2
- MMseqs2
- MDTraj 1.9.9
- scikit-learn 1.2.2
- ChimeraX 1.7

Python dependencies are listed in `requirements.txt`. A conda environment file is provided as `environment.yml`.

## Installation

Clone the repository:

    git clone https://github.com/JingyiYe/DNA-Origami-SMFS-CaE-analysis.git
    cd DNA-Origami-SMFS-CaE-analysis

Create the Python environment:

    conda env create -f environment.yml
    conda activate smfs-cae-analysis

Alternatively, install the Python dependencies using:

    pip install -r requirements.txt

This environment supports the Python-based HMM fitting interface and AF-Cluster post-processing scripts. MATLAB, LabVIEW, ChimeraX, ProteinTools, ColabFold and MMseqs2 should be installed or accessed separately as described in the protocol and documentation.

## Runnable example

A representative AF-Cluster example dataset, including CaE input files, AF-Cluster/ColabFold output files and execution instructions, is deposited in figshare:

**figshare DOI:** [10.6084/m9.figshare.32861552](https://doi.org/10.6084/m9.figshare.32861552)

The example is intended to demonstrate the analysis workflow and file formats. It does not include the complete raw optical-tweezers dataset.

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
