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

Raw and derived data are temporarily included in data/ for review.

Once the data paper is published, this repository will be updated with the DOI link.

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

