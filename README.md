# Cost Allocation in Optimized PyPSA Networks


This repository contains the workflow used to create figures and tables in **"Tracing Prices: a flow-based cost allocation for optimized networks"** ([see preprint](https://arxiv.org/abs/2010.13607))



## Initialize 


<!-- **For linux and mac:** clone the repository and run the install.sh script in your shell.  -->
For initializing clone the repo and create the python environment

```
git clone https://github.com/PyPSA/pypsa-eur.git ./pypsa-de
conda env create -f environment.yaml
```

## Run

To run the `snakmake` workflow make sure that the environment is activated 
```
conda activate costallocation
```
and start the run 
```
snakemake
```

### Further configuration 

Further instructions for modifying and extending will follow soon.


## License

This work is licensed under a [Creative Commons Zero v1.0 Universal](https://creativecommons.org/choose/zero/)
