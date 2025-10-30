# RangeX: Integral Projection Analyses


This repository contains cleaned and reproducible code accompanying the paper:

> **Dispersal limitation and competition explain lags in plant elevational range shifts**  
> Evelin Iseli, Anne Bjorkman, Vincent Ralph Clark, Anna Hargreaves, Paul Kardol, Vigdis Vandvik, Shengwei Zong, Jake Alexander  
> Submitted to _Science_, 2025

Code & repository by: Evelin Iseli

---

### Repository structure

- code/ : all analysis scripts (data prep, main analysis, cluster code, functions)
- data/ : raw and derived data used in analyses
- plots/ : figures and visualizations produced by the scripts
- output/ : result tables and model outputs

---

### Reproducibility

This project uses **renv** for dependency management.  

To reproduce the same R environment run:

```r
install.packages("renv")
renv::restore()
```

Then run the scripts in this order:
1. code/RangeX_IPM_dataprep_20251027.R
2. code/RangeX_IPM_maincode_20251027.R
3. code/RangeX_IPM_clustercode_20251027.R (optional --> the output of the cluster code is stored in data/derived and loaded in RangeX_IPM_maincode_20251027.R)

---

### Data availability

This repository includes all code necessary to reproduce the main IPM analyses. However, several large raw or intermediate datasets have been excluded from the GitHub repository due to size limits (>100 MB) or pending publication.

#### Included
- Cleaned demographic and germination datasets used for IPM fitting (`data/raw/`).
- Plots and IPM output summaries.

#### Excluded (too large for GitHub)
- High-frequency environmental log data (`RangeX_clean_EnvTMS4_2021_2023_CHE.csv`, `RangeX_clean_EnvHOBO_2021_2023_CHE.csv`).
- Downscaled macroclimatic temperature data based on CHELSAcruts and TerraClimate (downscaling as in Iseli et al. (2025), https://doi.org/10.1111/1365-2745.70114).
- Species distribution data for the study area (publicly available in aggregated form via InfoFlora (https://www.infoflora.ch)).
- Large bootstrapped cluster input and output files (`IPM_BootparaLong_mixed_*.csv`, `IPM_bootstrappedLambda_para_cluster_*.csv`, `IPM_bootstrappedVR_para_cluster_*.csv`).

These files are archived locally and will be deposited in a public data repository (e.g. Dryad or Zenodo) upon publication.  
Once available, this README will be updated with permanent DOI links and instructions for full reproducibility.

#### Reproducibility note
All intermediate IPM results can be reproduced from the included code, though some steps (e.g. model bootstrapping and full IPM fitting) are computationally intensive and were run on a cluster environment.

Note that one preprocessing script (`code/RangeX_dataprep_envdat_20251027.R`, part 2) requires local environmental and species distribution datasets that are not included in this repository due to data-sharing and size restrictions. This step prepares environmental summaries but is not required to reproduce the IPM analyses or population growth rate estimates. 

---

### Citation

If you use this code, please cite:

Evelin Iseli et al. (2025). RangeX Integral Projection Model Analyses.

GitHub repository: https://github.com/EvelinIseli/RangeX_IPMs_public

---

## Contact

For questions or collaborations, contact:
Evelin Iseli, evelin.iseli@bluewin.ch

---

## License

This repository is shared under the [CC BY-NC-ND 4.0 license](LICENSE).  
You are welcome to view and download the materials, but reuse, modification, or redistribution of the code or data is not permitted without permission.

