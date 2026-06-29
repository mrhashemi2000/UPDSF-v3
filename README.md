[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.20988680-blue)](https://doi.org/10.5281/zenodo.20988680) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)


# UPDSF-v3
The UPDSF v3.0 is a computational simulation engine designed to investigate the prebiotic transition from RNA-like nucleotides to DNA-like nucleotides.

## Author: Seyed Mohammad Reza Hashemi (Reza Hashemi)

Environment: 🐍 Python 3.8+

https://doi.org/10.5281/zenodo.20825578

https://doi.org/10.5281/zenodo.20733760 

https://doi.org/10.5281/zenodo.20759622 

https://doi.org/10.5281/zenodo.20771213 

https://doi.org/10.5281/zenodo.18594133

## Project Overview
The UPDSF v3.0 is a high-fidelity simulation engine designed to model the chemical selection of DNA nucleotides under prebiotic conditions. The framework specifically analyzes the preferential enrichment of Thymine (T) over Uracil (U) by calculating temperature-dependent (Arrhenius) and pH-dependent hydrolysis and polymerization rates.

## Scientific Basis
The simulation is calibrated based on established prebiotic chemistry literature:
- Ferris (1996): Polymerization kinetics.
- Huang (2005) & Joshi (2018): Clay surface catalysis and stability ratios.
- Cleaves (2010): pH-dependent protonation states and hydrolysis.

## Key Features
- 2D Sensitivity Analysis: Optimization of Temperature (60°C - 100°C) and pH (5.0 - 10.0).
- Kinetic Modeling: Implementation of Arrhenius equations with pH modulation.
- Clay Protection Factor (CPF): Modeling the protective effect of mineral surfaces on DNA stability.
- Comprehensive Visualization: Generation of 3D response su

## Installation
Clone the repository and install the dependencies:
git clone https://github.com/mrhashemi2000/UPDSF-v3.git
cd UPDSF-v3

pip install -r requirements.txt

## Usage
Run the main simulation engine to find the optimal prebiotic conditions:
python src/main.py

## Citation 

If you use this framework in your research, please cite it as: Hashemi, S. M. R. (2026). Unified Prebiotic DNA Selection Framework (UPDSF) v3.0 https://doi.org/10.5281/zenodo.20988680

## References
The model is calibrated based on the following landmark studies in prebiotic chemistry:

Ferris, J. P., Hill, A. R., Liu, R., & Orgel, L. E. (1996)
Synthesis of long prebiotic oligomers on mineral surfaces.
Nature, 381(6577), 59–61.
DOI: 10.1038/381059a0

Cleaves, H. J. (2010)
The origin of the biologically coded amino acids.
Journal of Theoretical Biology, 263(4), 490–498.
DOI: 10.1016/j.jtbi.2009.12.014
