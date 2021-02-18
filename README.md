## Spline-based accelerated failure time model

We provide a few scripts to estimate spline-based accelerated failure time 
models as suggested by Pang et al. (2021). 
The examples can be run by running `R CMD BATCH --no-restore --no-save run.R` 
which produces the output `run.Rout` and the plots. The `fit-saft.R` file loads 
the packages that are needed and assigns a function to estimate the model. 

The estimation method does not provide gradients and therefore the estimation 
time might be greatly reduced using the gradients as well. 

## References
Pang, M, Platt, RW, Schuster, T, Abrahamowicz, M. 2021. “Spline‐based 
accelerated failure time model.” *Statistics in Medicine* 40: 481– 497.
<https://doi.org/10.1002/sim.8786>.
