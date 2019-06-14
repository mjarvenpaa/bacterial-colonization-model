# Readme

This readme contains a brief introduction to the code used in the simulation experiments and analysis of MRSA data in the paper *A Bayesian model of acquisition and clearance of bacterial colonization incorporating within-host variation* by Järvenpää et al. (2019) PLoS Computational Biology, 15(4):e1006534

Link to the paper: <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006534>. 

## List of code files in /code/:

- `mixture_model_demos.R`: Scripts for running the simulation experiments with the mixture model in the paper.


- `mixture_model_diags.R`: Functions for assessing the convergence of MCMC, computing final estimates from the MCMC output, posterior predictive checks and related visualisations of the fitted model.

- `fit_mixture_model.R`: The Gibbs sampling algorithm for fitting the mixture model.

- `mixture_model.R`: Various help functions and visualisations for mixture model fitting. 


- `simulation_functions.R`: Functions for running the Wright-Fisher simulation model.  


- `cond_dens_of_sim_model_for_gibbs`: Function for running the simulations that are used to estimate the conditional densities p\_S(d\_1i|n\_eff,mu).

- `cond_dens_help_functions.R`: Help functions needed for the computation of the conditional densities p\_S(d\_1i|n\_eff,mu) and some visualisations.


- `data_summary_functions.R`: Functions for computing distances, summaries and the discrepancy values needed for the ABC inference.

- `abc_simul_cluster.R`: Function for running the simulations and computing the simulated discrepancy functions for ABC inference.

- `abc_methods_simul_model.R`: Functions for ABC-inference (and related visualisations).


- `figure_functions.R`: Misc. plotting functions.

- `auxiliary_functions.R`: Misc. help function.

- `get_settings.R`: General settings and parameter values etc. for the data analysis are gathered here.

- `get_data.R`: Returns the data sets.


## List of data files:

We provide the same data D used for mixture model fitting in the article. Note that data set contains only the distance and time difference pairs necessary for fitting the mixture model. The raw data is not provided here.

- `/data_CLEAR_minimal/mrsa_clear_testdata_armE.RData`: Data from the 'education' (E) arm of the CLEAR study. 

We do not provide the external data D\_0 used for ABC analysis here but, instead, we provide similar 'artificial' data sets generated from the model itself:

- `/data_external/distance_test_data.RData`: Test data set that can be used for ABC analysis. 


## Precomputed simulations:

Simulations needed for ABC analysis and for computing the conditional densities are costly and were thus precomputed in a computer cluster system. Luckily, fitting the mixture model is relatively fast because these costly precomputed simulations can be re-used when fitting the mixture model to different data sets. We provide these precomputed values in the following R data files:

- `/simulation_outputs_abc/abc_final_simul/abc_results/ABC_samples_l1.RData`: Contains the ABC posterior of \mu,n\_eff given data D_0 computed in 50x50 grid of parameter values. 

- `/simulation_outputs_cond_dens/cond_dens_final_simul/precomp_density.RData`: Contains the precomputed estimates of conditional densities in 100x100 grid of parameter values.


## Additional R packages required:

- latex2exp (for plotting)
- coda (for assessing the convergence of MCMC)
- BayesianTools (for plotting)
- MASS (for kernel density estimation used in the plotting)

These packages can be installed in a usual way i.e. from R terminal using `install.packages()` command.


## Instructions to get started:

Once the code files are placed to a desired location and the required extra R packages installed, the mixture model can be fitted to the provided MRSA data using `compute.true.data` function and to 'artificial data' using `compute.demos` function in the file `mixture_model_demos.R`. Note that these computations take some time and output files together require ~1Gb of hard disk space. The user needs to provide the filepath to the main folder (i.e. the folder containing 'code' subfolder) as the input variable 'root'. The output files and visualisations are saved to the subfolder `simulation_outputs_mixture_model`. 

Specifically:

- `compute.true.data`: This script fits the mixture model to the MRSA data set and generates the Figures 4,7 and 8 in the article (as well as some extra figures).

- `compute.demos`: This script fits the mixture model to multiple 'artificial' data sets that are generated from the model itself to assess the correctness and efficiency of the parameter fitting algorithm. Figure 6 of the paper is generated as output. 

Note:
We provide easy-to-use code for mixture model fitting (check out the two function above). Full scripts for precomputing the conditional densities and running the ABC analysis which are necessary for fitting the mixture model are also provided but these have high computational cost. To replicate these analyses one may need to modify these code files so that these can be run on a computer cluster. The precomputed values resulting from such computation are provided as R data files as described above. Note also that the data D\_0 could not be provided here but we instead provide similar 'fake' data set. 


## Support:

If you have questions about the code or the related experiments, please contact marko.j.jarvenpaa<at>aalto.fi 


## License:

This code is provided under the MIT licence.


