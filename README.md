# runoutSIM
The `runoutSIM` R package is developed for regional-scale runout simulations using random walks. It currently can be used to estimate the spatial extent, velocity (via [Perla et al's](https://www.cambridge.org/core/journals/journal-of-glaciology/article/twoparameter-model-of-snowavalanche-motion/B87923FFC6ADAF61B0079EEBCBD96F19) two-parameter friction model), and connectivity of:
* Debris flows
* Snow avalanches

**Features**
* Random walk simulation with slope-based transition probabilities
* Velocity modeling using Perla et al.'s (1980) friction law
* Connectivity analysis to assess impact on downslope features
* Allows easy workflow for source area prediction models (e.g., GAM, machine learning)
* Interactive mapping with leaflet and htmlwidgets
* Optimized for parallel processing support for large-scale simulations
* Integrated into ['runoptGPP'](https://github.com/jngtz/runoptGPP) for grid search parameter optimization


**Installation**

You can install runoutSIM with:

```r
devtools::install_github("jngtz/runoutSIM")
```
or 

```r
remotes::install_github("jngtz/runoutSIM")
```

## Contributing
We welcome contributions! Please open an issue or submit a pull request. For major changes, start by discussing your ideas via an issue.

## Citation
If you use runoutSIM in your research, please cite: 

Goetz, J. (2025). runoutSIM: An open-source R package for simulating mass movement runout and connectivity using random walks. GitHub Repository

