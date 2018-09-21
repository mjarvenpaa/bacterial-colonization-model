# This function can be used for running the (computationally costly) simulations and computing the 
# discrepancies for the ABC analysis of the data 'D_0'. The idea is that this function is first used to
# compute all the required discrepancies in a computer cluster and the actual ABC posterior for 
# parameters \mu and n_eff, which is used as a prior density in mixture model, as well as some other 
# analysis of the simulations, can be then computed easily on a single computer.
#
# Note: Computation time is very high. Computing is however embarrassingly parallelizable 
# and full simulation is here already splitted to smaller tasks so that each input 'id' 
# (integer from 1 to 1000) represents one computation task. That is, to compute the full 
# set of simulations one needs to run the function below in a loop that goes through all
# 'id' values from 1 to 1000. 
# 
# Note that the computational cost of running these simulations is very high so to estimate 
# these densities, one should run the code below in a computer cluster environment. 

run.abc.simul <- function(root, id) {
  
  simulation.name <- 'abc_final_simul' 
  
  # NOTE: 'simulation+distances' only saves the distance matrices to file
  # This is useful because full simulation outputs can require huge amount of hard disk space!
  methods.to.run <- c('simulation+distances','plot.distance.evolution','discrepancies')
  
  always.recompute <- TRUE # if TRUE then the results are recomputed even if corresponding file already exists
  
  ###########################################
  
  ## libraries, functions
  code.path <- paste(root, '/code', sep='')
  source(paste(code.path, '/get_settings.R', sep=''))
  source(paste(code.path, '/get_data.R', sep=''))
  source(paste(code.path, '/simulation_functions.R', sep=''))
  source(paste(code.path, '/data_summary_functions.R', sep=''))
  source(paste(code.path, '/abc_methods_simul_model.R', sep=''))
  source(paste(code.path, '/auxiliary_functions.R', sep=''))
  source(paste(code.path, '/figure_functions.R', sep=''))
  
  # Create a directory for the simulation if it does not yet exist:
  res.dir.name <- paste(root, '/simulation_outputs_abc/', simulation.name, sep='')
  if (!file.exists(res.dir.name)) {
    dir.create(res.dir.name)
    dir.create(paste(res.dir.name, '/saved_populations', sep=''))
    dir.create(paste(res.dir.name, '/distances', sep=''))
    dir.create(paste(res.dir.name, '/discrepancies', sep=''))
    dir.create(paste(res.dir.name, '/figures', sep=''))
    dir.create(paste(res.dir.name, '/abc_results', sep=''))
  }
  
  ###########################################
  
  ## Get the settings for the current run
  opt <- get.abc.fit.settings(simulation.name=simulation.name)
  
  # number of parallel runs and maximum command line argument
  # THIS COULD BE SMALLER OR THE SAME AS N.REP WHICH IS THE TOTAL # OF REPETITIONS FOR EACH PARAMETER
  parallel.runs <- opt$parallel.runs
  
  ## based on the number of parallel runs, choose the repetition indexes for this run
  n.rep <- opt$n.rep
  reps <- get.current.ids(id, parallel.runs, n.rep)
  print.output <- opt$print.output
  
  genome.length <- opt$genome.length
  
  n1 <- opt$n1
  n2 <- opt$n2
  neff.grid <- opt$neff.grid
  mu.grid <- opt$mu.grid
  
  
  ## load true observed data to get the number of repeated simulations to run and the corresponding time points
  obs.data <- get.sim.model.data(root=root)
  data <- list(first=obs.data$data.first, last=obs.data$data.last)
  data.AM <- obs.data$data.AM
  
  
  # number of simulations for matching the observed data i.e. the number of patients (the 8 (or 7) patients + patients A-M (13))
  nr.data.sim <- length(data$first)
  nr.data.sim.AM <- length(data.AM)
  
  save.interval <- opt$save.interval
  save.interval.AM <- opt$save.interval.AM
  
  # For the discrepancy computation, which discrepancy/discrepancies to compute
  discrepancies <- opt$discrepancies
  
  
  ########################
  # Run simulation model #
  ########################
  
  pops.filename.id <- paste(res.dir.name, '/saved_populations/pops_id', id, '.RData', sep='')
  dist.filename.id <- paste(res.dir.name, '/distances/dists_id', id, '.RData', sep='')
  save.pops <- !any(methods.to.run == 'simulation+distances') # whether pops are saved
  
  if (any(methods.to.run == 'simulation') || any(methods.to.run == 'simulation+distances')) {
    
    simul.save.freq <- 1
    dist.save.freq <- 1
    
    # check first if some computations already exist
    if (always.recompute || !file.exists(pops.filename.id)) {
      if (!save.pops) {
        dist.data <- array(list(), dim=c(n1,n2,length(reps),nr.data.sim))
        dist.data.AM <- array(list(), dim=c(n1,n2,length(reps),nr.data.sim.AM))
      } else {
        sim.data <- array(list(), dim=c(n1,n2,length(reps),nr.data.sim)) # index (mu,n_eff,rep,nr) and data is pops format which is a list
        sim.data.AM <- array(list(), dim=c(n1,n2,length(reps),nr.data.sim.AM))
      }
      
      for (i in 1:length(reps)) {
        # seed is set according to the current rep index
        set.seed(123+reps[i])
        
        for (j in 1:length(neff.grid)) { 
          for (k in 1:length(mu.grid)) { 
            
            ## simulate with current (neff,mu)-parameter
            if (prior.neff.mu.value(neff.grid[j], mu.grid[k], neff.grid, mu.grid) == 0) {
              if (print.output) {
                print(paste('Prior density is zero, neff=',neff.grid[j],', mu=',mu.grid[k],sep = ''))
              }
              next
            }
            if (print.output) {
              print(paste('Simulating with neff=',neff.grid[j],', mu=',mu.grid[k],sep = ''))
            }
            
            ## DATA 1-8
            pops <- rep(list(),nr.data.sim)
            for (l in 1:nr.data.sim) {
              if (print.output) { print(paste('l=',l,'/',nr.data.sim,sep = '')) }
              
              # returns pops which contains all the generations for current (mu,n_eff) and rep
              pops[[l]] <- simulate.evolution(mutation.rate=mu.grid[k], n.strains=neff.grid[j], genome.length=genome.length, 
                                              n.generations.to.simulate=save.interval[length(save.interval)], save.interval=save.interval, 
                                              save.path=save.path, save.extension=save.extension, print.output=print.output, save.output=FALSE)
            }
            # add new output to the datastructure
            if (save.pops) {
              sim.data[j,k,i,] <- pops
            }
            
            ## DATA A-M
            pops.AM <- rep(list(),nr.data.sim.AM)
            for (l in 1:nr.data.sim.AM) {
              if (print.output) { print(paste('l=',l,'/',nr.data.sim.AM,sep = '')) }
              
              # as above, but patients A-M
              pops.AM[[l]] <- simulate.evolution(mutation.rate=mu.grid[k], n.strains=neff.grid[j], genome.length=genome.length, 
                                                 n.generations.to.simulate=save.interval.AM[length(save.interval.AM)], save.interval=save.interval.AM, 
                                                 save.path=save.path, save.extension=save.extension, print.output=print.output, save.output=FALSE)
            }
            # add new output to the datastructure
            if (save.pops) {
              sim.data.AM[j,k,i,] <- pops.AM
            }
            
            # Compute also distances now here!
            ##################################
            if (!save.pops) {
              for (l in 1:nr.data.sim) {
                # resampling to the same size as the observed data
                # NOTE: We assume here that the first and second time point have the same size which holds here
                n.strains.new <- get.nr.individuals.from.nr.pairs(length(data$first[[l]]))
                pops.l <- resample.pop(pops=pops[[l]], n.strains.new=n.strains.new)
                
                # compute the distances for each generated population 
                for (m in 1:length(pops.l)) {
                  dist.data[[j,k,i,l]][[m]] <- compute.distance.distribution(pops.l[[m]])
                }
              }
              for (l in 1:nr.data.sim.AM) {
                # as above, but patients A-M
                n.strains.new <- get.nr.individuals.from.nr.pairs(length(data.AM[[l]]))
                pops.AM.l <- resample.pop(pops=pops.AM[[l]], n.strains.new=n.strains.new)
                
                for (m in 1:length(pops.AM.l)) {
                  dist.data.AM[[j,k,i,l]][[m]] <- compute.distance.distribution(pops.AM.l[[m]])
                }
              }
            }
          }
        }
        
        # save populations in the datastructures to file every now and then (overwriting the existing one)
        if (save.pops && (i == 1 || i %% simul.save.freq == 0 || i == length(reps))) {
          save(sim.data, sim.data.AM, file=pops.filename.id)
        }
        
        # save distance if they were computed already here
        if (!save.pops && (i == 1 || i %% dist.save.freq == 0 || i == length(reps))) {
          save(dist.data, dist.data.AM, file=dist.filename.id)
        }
      }
    } else if (print.output) {
      print('Simulation datafile already exists.')
    }
  }
  
  #####################################################
  # Compute distance distributions at given intervals #
  #####################################################
  
  if (save.pops && any(methods.to.run == 'distances')) {
    
    dist.save.freq <- 1
    
    # check first if some computations already exist
    if (always.recompute || !file.exists(dist.filename.id)) {
      
      dist.data <- array(list(), dim=c(n1,n2,length(reps),nr.data.sim))
      dist.data.AM <- array(list(), dim=c(n1,n2,length(reps),nr.data.sim.AM))
      
      # load pops data (unless it was just computed)
      if (!exists('sim.data') || !exists('sim.data.AM')) {
        load(file=pops.filename.id) # sim.data, sim.data.AM
      }
      
      # compute the distances
      for (i in 1:length(reps)) {
        # seed is set according to the current rep index
        set.seed(123+reps[i])
        
        for (j in 1:length(neff.grid)) { 
          for (k in 1:length(mu.grid)) {
            if (prior.neff.mu.value(neff.grid[j], mu.grid[k], neff.grid, mu.grid) == 0) {
              next
            }
            
            for (l in 1:nr.data.sim) {
              # resampling to the same size as the observed data
              # NOTE: We assume here that the first and second time point have the same size which holds here
              n.strains.new <- get.nr.individuals.from.nr.pairs(length(data$first[[l]]))
              pops <- resample.pop(pops=sim.data[[j,k,i,l]], n.strains.new=n.strains.new)
              
              # compute the distances for each generated population 
              for (m in 1:length(pops)) {
                dist.data[[j,k,i,l]][[m]] <- compute.distance.distribution(pops[[m]])
              }
            }
            for (l in 1:nr.data.sim.AM) {
              # as above, the second data set 
              
              n.strains.new <- get.nr.individuals.from.nr.pairs(length(data.AM[[l]]))
              pops.AM <- resample.pop(pops=sim.data.AM[[j,k,i,l]], n.strains.new=n.strains.new)
              
              for (m in 1:length(pops.AM)) {
                dist.data.AM[[j,k,i,l]][[m]] <- compute.distance.distribution(pops.AM[[m]])
              }
            }
          }
        }
        
        # save the distances to file every now and then (overwriting the existing one)
        if (i == 1 || i %% dist.save.freq == 0 || i == length(reps)) {
          save(dist.data, dist.data.AM, file=dist.filename.id)
        }
      }
    } else if (print.output) {
      print('Distance datafile already exists.')
    }
  }
  
  
  #######################################################################
  # Plot distance evolution, this is only for illustrations and testing #
  #######################################################################
  if (any(methods.to.run == 'plot.distance.evolution') && id == 1) {
    # plots now only for the first repetition and (!) in the first batch run 
    
    # generations to be plotted 
    if (length(save.interval) == 1) {
      generations.to.consider.plotting <- seq(1,n.generations.to.simulate,by=save.interval)
    } else {
      generations.to.consider.plotting <- save.interval
    }
    figure.path <- paste(res.dir.name, '/figures', sep='')
    
    # load distance data
    load(file=dist.filename.id) # dist.data, dist.data.AM
    
    for (i in c(1)) { # plot only first repetition simulation in the parallel task for now
      for (j in 1:length(neff.grid)) { 
        for (k in 1:length(mu.grid)) {
          if (prior.neff.mu.value(neff.grid[j], mu.grid[k], neff.grid, mu.grid) == 0) {
            next
          }
          
          for (l in c(1)) { # plot only the simulation corresponding the first observed data point for the first data set for now
            
            save.extension <- formulate.save.extension(c('neff'=neff.grid[j], 'mu'=mu.grid[k], 'rep'=reps[i], 'obs_sim'=1))
            dists <- dist.data[[j,k,i,l]]
            plot.distance.evolution.pops(dists, generations.to.consider.plotting, save.extension, figure.path, adapt.bins = TRUE)
          }
        }
      }
    }
  }
  
  
  ##################################################################
  # Compute discrepancies for each replication/simulated parameter #
  ##################################################################
  if (any(methods.to.run == 'discrepancies')) {
    
    # Computes now all the provided discrepancies
    # One can then compare the results with different dicrepancies without recomputing stuff.
    for (discrepancy.ind in 1:length(discrepancies)) {
      
      discrepancy.name <- discrepancies[[discrepancy.ind]]$name
      cur.discr.filename.id <- paste(res.dir.name, '/discrepancies/', discrepancy.name, '_id', id, '.RData', sep='')
      
      # load the distance data
      discrepancy.data <- array(NA,dim=c(n1,n2,length(reps)))
      load(file=dist.filename.id)
      
      # compute the discrepancy for each simulation (computing this is fast)
      for (i in 1:length(reps)) {
        # seed is set according to the current rep index
        set.seed(123+reps[i])
        
        for (j in 1:length(neff.grid)) { 
          for (k in 1:length(mu.grid)) {
            
            if (prior.neff.mu.value(neff.grid[j], mu.grid[k], neff.grid, mu.grid) == 0) {
              # set discrepancy always to inf if outside prior support
              discrepancy.data[j,k,i] <- Inf
              next
            }
            dists <- dist.data[j,k,i,]
            dists.AM <- dist.data.AM[j,k,i,]
            discrepancy.data[j,k,i] <- discrepancy(discrepancies[[discrepancy.ind]], dists, dists.AM, data, data.AM, save.interval, save.interval.AM, 
                                                   print.discrepancy=print.output, neff=neff.grid[j], mu=mu.grid[k]) 
          }
        }
      }
      
      # save the discrepancies
      save(discrepancy.data, file=cur.discr.filename.id)
      
      # debug printing for discrepancies
      if (print.output) {
        print(paste(discrepancy.name, ':', sep = ''))
        print(discrepancy.data)
      }
    }
  }
}




