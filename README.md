# ABL resistance
This repository will be used to benchmark and improve KinoML and Perses.

### Ideas

- reproduce ABL:inhibitor structures from [Hauser 2018](https://doi.org/10.1038/s42003-018-0075-x) using the KinoML pipeline
- run mutation benchmark from [Aldeghi 2019](https://pubs.acs.org/doi/abs/10.1021/acscentsci.9b00590) using OpenFF 1.2.0
- absolute free energy calculations with Yank might be interesting too
- long simulations to analyze stability of a variety of ABL:inhibitor complexes
- dock ATP/Mg2+ in binding pocket and analyze the effect of point mutations
- scale up to KINOMEScan data

### How to use this repository

1. Clone repository

`git clone https://github.com/openkinome/abl_resistance`

2. Create Conda environment
  
`conda env create -f environment.yml`  
`conda activate abl_resistance`

### Structure

- `notebooks/atp_kinase_conformations.ipynb`  
  - jupyter notebook analyzing the conformations of ATP bound kinases
- `notebooks/abl1_atp_modeling.ipynb`
  - jupyter notebook generating an ABL1 ATP complex
- `notebooks/abl_complex_modeling.ipynb`
  - jupyter notebook generating inhibitor bound complexes for the [Hauser 2018](https://doi.org/10.1038/s42003-018-0075-x) benchmark
- `notebooks/abl_complex_modeling_with_water.ipynb`
  - jupyter notebook using updated KinoML functionalities to generate inhibitor bound complexes for the [Hauser 2018](https://doi.org/10.1038/s42003-018-0075-x) benchmark including water

### Authors

- David Schaller <david.schaller@charite.de>
- William Glass <william.glass@choderalab.org>

### License
This repository is licensed under the MIT license.