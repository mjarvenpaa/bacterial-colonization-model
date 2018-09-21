# This function can be used for computing the estimates of the intractable conditional 
# densities p(d_i1|mu,n_eff) by raw simulations. These estimated densities are needed
# when fitting the mixture model.
#
# Note: Computation time is very high. Computing is however embarrassingly parallelizable 
# and full simulation is here already splitted to smaller tasks so that each input 'id' 
# (integer from 1 to 1000) represents one computation task. That is, to compute the full 
# set of simulations one needs to run the function below in a loop that goes through all
# 'id' values from 1 to 1000. 
# 
# Note that the computational cost of running these simulations is very high so to estimate 
# these densities, one should run the code below in a computer cluster environment. 

run.cond.dens.simul <- function(root, id) {

  simulation.name <- 'cond_dens_final_simul'
  
  ## libraries, functions
  code.path <- paste(root, '/code', sep='')
  source(paste(code.path, '/get_settings.R', sep=''))
  source(paste(code.path, '/abc_methods_simul_model.R', sep=''))
  source(paste(code.path, '/data_summary_functions.R', sep=''))
  source(paste(code.path, '/simulation_functions.R', sep=''))
  source(paste(code.path, '/auxiliary_functions.R', sep=''))
  
  ####################################################################
  # Compute the densities p(d_1i|n_eff,mu) for mixture model fitting #
  ####################################################################
  
  ## Go through the different mu values in the mu.grid
  d1i.save.freq <- 5 # how often the results are saved to file
  
  opt2 <- get.cond.dens.settings(simulation.name=simulation.name)
  
  print.output <- opt2$print.output
  
  # these are now saved to different location as original simulations
  cond.dens.dir.name <- paste(root, '/simulation_outputs_cond_dens/', simulation.name, sep='')
  cond.dens.filename <- paste(cond.dens.dir.name, '/cond_dens_id', id, '.RData', sep='')
  
  if (!file.exists(cond.dens.dir.name)) {
    dir.create(cond.dens.dir.name)
  }
  
  genome.length <- opt2$genome.length
  
  parallel.runs <- opt2$parallel.runs
  nr.samples <- opt2$nr.samples 
  save.interval <- opt2$save.interval 
  simul.save.interval <- opt2$simul.save.interval
  n.generations.to.simulate <- save.interval[length(save.interval)]
  
  
  # grid for cond dens simulations
  neff.grid.new <- opt2$neff.grid 
  mu.grid.new <- opt2$mu.grid 
  
  # NOTE: save.interval is when the populations are computed and returned, simul.save.interval is
  # when the distances are computed and saved.
  # In practice these two are set to be the same (or at least simuls.save.interval
  # must be subset of save.interval)
  
  ## based on the number of parallel runs, choose the repetition indexes for this run
  reps <- get.current.ids(id, parallel.runs, nr.samples)
  
  # initialize
  d_i1.samples <- array(NA, c(length(neff.grid.new), length(mu.grid.new), length(reps), length(simul.save.interval)))
  
  for (i in 1:length(reps)) {
    # seed is set according to the current rep index
    set.seed(1234+reps[i])
    
    for (j in 1:length(neff.grid.new)) {
      for (k in 1:length(mu.grid.new)) {
        
        # check prior bound
        if (prior.neff.mu.value(neff.grid.new[j], mu.grid.new[k], neff.grid.new, mu.grid.new) == 0) {
          if (print.output) {
            print(paste('Prior density is zero, n_eff=',neff.grid.new[j],', mu=',mu.grid.new[k],sep = ''))
          }
          next
        }
        
        if (print.output) {
          print(paste('Simulating with n_eff=', neff.grid.new[j], ', mu=', mu.grid.new[k], sep = ''))
        }
        # Simulate new d_i1 data from the simulation model with the given (n_eff,mu)-parameter
        pop.j <- simulate.evolution(mutation.rate=mu.grid.new[k], n.strains=neff.grid.new[j], genome.length=genome.length, 
                                    n.generations.to.simulate=n.generations.to.simulate, save.interval=save.interval, 
                                    save.path=save.path, save.extension=save.extension, print.output=print.output, save.output=FALSE)
        
        for (l in 1:length(simul.save.interval)) {
          # Pick the d_i1 value from the generated population randomly
          # Select the two genomes randomly (at the same time point) and compute their distance
          # TODO: could this be repeated also for each sampled neff,mu -> less simulations needed!
          save.ind <- which(simul.save.interval[l]==save.interval)
          if (length(save.ind) != 1) {
            stop('Something is wrong with pop indexes.')
          }
          d_i1.samples[j,k,i,l] <- random.pop.distance(pop = pop.j[[save.ind]], n = 1)
        }
        
        # For each n_eff,mu in grid we have now the density p(d_1i|n_eff,mu) estimated by sampling.
        # This density is represented by non-negative integer-valued samples d_ij for each grid point ((n_eff)_j,mu_j)
        # and, possibly, at some different simul. time points
      }
    }
    
    ## Save estimated densities to file (to be later used for the Gibbs sampling algorithm)
    if (i == 1 || i %% d1i.save.freq == 0 || i == length(reps)) {
      # It is convenient to save also the grid
      save(d_i1.samples, neff.grid.new, mu.grid.new, file = cond.dens.filename)
    }
  }
}




