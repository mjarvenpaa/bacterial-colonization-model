# Contains functions that return settings for different experiments.
#
# All the important settings etc. are gathered here so that the experiments can be to rerun easily.
# If one needs to know some used setting later on, one can find it easily here.
#
# The basic idea of the settings is the following: The basic inference phases are: 
# 1a) [cluster] run simulations in a (neff,mu)-parameter grid and compute the (ABC) discrepancies using data D_0
# 1b) compute the ABC posterior p_epsilon(mu,n_eff|D_0) using the above model simulations
# 2a) [cluster] use extra simulations in the same or different grid as 1) to estimate densities p(d_1i|mu,n_eff)
# 2b) estimate mixture model parameters 'theta' using the output of 1) and 2a) and additional data D i.e. CLEAR data
#
# Each of these tasks features some settings (decision choices etc.) and depends on the previous
# computations and thus also settings used to do those previous computations.

## Settings for phase 1:

get.abc.fit.settings <- function(simulation.name, print.settings=TRUE) {
  
  opt <- list() # settings are saved to this list
  
  opt$simulation.name <- simulation.name
  
  opt$parallel.runs <- 1
  opt$print.output <- TRUE
  
  ## General settings
  ###################
  
  opt$genome.length <- 3000000
  
  # when to compute&save the summaries etc. (not all of these values are needed for ABC analysis)
  opt$save.interval <- unique(c(1, seq(200,1000,by=200), seq(1000,6000,by=500), seq(6000,9000,by=500)))
  opt$save.interval.AM <- unique(c(1, seq(500,9000,by=500)))
  
  time1 <- 1000 # ~ 2 months, we actually only know that this is an upper bound
  time2 <- 6000 # ~ 1 year
  time.AM <- 6000 # ~ 1 year, these are actually not known
  
  times1 <- rep(time1, 8) # we assume same measurement times for each data points in 1-8 and also for A-M
  times2 <- rep(time2, 8) 
  times.AM <- rep(time.AM, 13)
  
  ## Which discrepancy/discrepancies to compute and related settings
  # One of these is finally chosen for further computations.
  # ignore1219=T if the patient 1219 is not used as data
  # ignoreAM=T if use only patients 1-8 as data (considered for testing only)
  opt$discrepancies <- list()
  opt$discrepancies[[1]] <- list(name='eucl_noAM',type='eucl', ignore1219=T,ignoreAM=T, times18_1=times1,times18_2=times2,timesAM=times.AM)
  opt$discrepancies[[2]] <- list(name='eucl',type='eucl', ignore1219=T,ignoreAM=F, times18_1=times1,times18_2=times2,timesAM=times.AM)
  opt$discrepancies[[3]] <- list(name='l1_noAM',type='l1_sep', ignore1219=T,ignoreAM=T, times18_1=times1,times18_2=times2,timesAM=times.AM)
  opt$discrepancies[[4]] <- list(name='l1',type='l1_sep', ignore1219=T,ignoreAM=F, times18_1=times1,times18_2=times2,timesAM=times.AM)
  
  ## ABC related settings: choosing the threshold
  opt$epsilon.strategy <- 'nr.samples<eps' # 'fixed', 'nr.samples<eps' or 'quantile'
  opt$epsilon.quantile <- 0.01
  opt$abc.nr.samples <- 5000
  opt$epsilon.fixed <- 4

  
  if (simulation.name == 'abc_test') {
    # for testing
    
    opt$parallel.runs <- 20
    opt$print.output <- FALSE
    opt$n.rep <- 20 # How many repeated simulations with each grid point (n_eff,mu)
    
    opt <- gen.neff.mu_grid(opt, grid.level='final20x20')
    
  } else if (simulation.name == 'abc_final_simul') {
    # FINAL RUN IN CLUSTER, we use 50x50 grid
    
    opt$parallel.runs <- 1000
    opt$print.output <- FALSE
    opt$n.rep <- 1000 # How many repeated simulations with each grid point (n_eff,mu)
    
    opt <- gen.neff.mu_grid(opt, grid.level='paper50x50')
    
  } else {
    stop(paste('Incorrect simulation.name in task 1: ', simulation.name, sep = ''))
  }
  
  if (print.settings) { print(opt) }
  return(opt)
}

####################################################################################
####################################################################################
## Settings for phase 2a:

get.cond.dens.settings <- function(simulation.name, print.settings=TRUE) {
  # Get settings for computing the conditional densities p(d_i1|mu,n_eff)
  
  ## settings saved to this list
  opt2 <- list()
  opt2$simulation.name <- simulation.name
  
  opt2$parallel.runs <- 1
  opt2$print.output <- TRUE
  
  opt2$genome.length <- 3000000
  
  # when the pairwise distances are computed (if considered fixed)
  # we actually may need only the last timepoint
  opt2$save.interval <- c(1000,2000,3000,4000,5000,6000)
  opt2$simul.save.interval <- opt2$save.interval # p(d_i1|mu,n_eff) is computed by simulation at these times 
  

  if (simulation.name == 'cond_dens_test') {
    # For testing
    
    opt2$parallel.runs <- 50
    opt2$print.output <- FALSE
    
    opt2$nr.samples <- 500 # How many simulations for each (mu,n_eff) in the grid
    opt2 <- gen.neff.mu_grid(opt2, grid.level='final20x20')
    
  } else if (simulation.name == 'cond_dens_final_simul') {
    ## We use finally 100x100 grid to avoid some artifacts due to not dense enough grid
    
    opt2$parallel.runs <- 1000
    opt2$print.output <- FALSE
    
    opt2$nr.samples <- 10000 # How many simulations for each (mu,n_eff) in the grid
    opt2 <- gen.neff.mu_grid(opt2, grid.level='paper100x100')
    
  } else {
    stop(paste('Incorrect simulation.name in task 2: ', simulation.name, sep = ''))
  }
  
  if (print.settings) { print(opt2) }
  return(opt2) 
}

####################################################################################
####################################################################################
## Settings for phase 2b:

get.mixture.fitting.settings <- function(simulation.name, print.settings=TRUE) {
  
  opt3 <- list()
  
  ## Which Gamma model, Exp if kt == 1 (this is the same as 'k' in the paper)
  opt3$kt <- 5
  
  ## Gibbs sampling settings
  opt3$nr.chains <- 3
  opt3$nr.gibbs.samples <- 20000
  opt3$rand.gibbs.init <- TRUE # whether to initialize by drawing from the prior
  
  #opt3$proposal.v.neff <- 200^2 # proposal variances in the Metropolis step
  #opt3$proposal.v.mu <- 0.00015^2
  opt3$proposal.v.neff <- 400^2 
  opt3$proposal.v.mu <- 0.00025^2
  
  
  ## Which ABC inference result is used in fitting the mixture model
  opt3$use.prior.from.abc <- 0
  opt3$chosen.discrepancy.name <- 'l1' # 'eucl_noAM', 'eucl', 'l1_noAM', 'l1'
  
  
  ## Hyperpriors
  opt3$alpha_lambda <- cg.params(10000, 10000^2, opt3$kt)[1]
  opt3$beta_lambda <- cg.params(10000, 10000^2, opt3$kt)[2]
  
  #opt3$gamma_omega <- c(0.5,0.5) # U-shaped
  opt3$gamma_omega <- c(1,1) # Uniform distribution
  
  ## Other settings
  opt3$print.output <- TRUE
  true.params <- NA
  
  
  if (simulation.name == 'mixture_test_artificial') {
    # For misc. testing
    
    # select the true parameters etc. for generating the 'fake' data
    true.params <- list()
    true.params$N <- 150
    true.params$true.neff.ind <- 22 # 22 is ~2000 in 100x100 grid
    true.params$true.mu.ind <- 22 # 22 is ~0.001 in 100x100 grid
    true.params$omega1 <- 0.7
    true.params$lambda <- 0.0002
    true.params$kt <- 5 # This is actually not estimated but assumed as fixed and known 
    true.params$hosp.times <- 100*c(5,10,15)
    
    opt3$use.prior.from.abc <- 0
    
    opt3$kt <- 5
    
    ## Gibbs sampling settings
    opt3$nr.chains <- 3
    opt3$nr.gibbs.samples <- 10000
    opt3$rand.gibbs.init <- TRUE # whether to initialize by drawing from the prior
    
    opt3$proposal.v.neff <- 200^2 # proposal variances in the Metropolis step
    opt3$proposal.v.mu <- 0.00015^2
    
    opt3$which.data <- 'artificial1'
    opt3$which.arm <- NA
    opt3$abc.simulation.name <- 'abc_final_simul'
    opt3$cond.dens.simulation.name <- 'cond_dens_final_simul'
    
    opt3$alpha_lambda <- cg.params(10000, 10000^2, opt3$kt)[1]
    opt3$beta_lambda <- cg.params(10000, 10000^2, opt3$kt)[2]
    
    
  } else if (substr(simulation.name,1,15) == 'artificial_demo') {
    # Format: artificial_demo_N=xxx_w=yyy_rep=zzz (where yyy is in procent e.g. yyy==15 corresponds omega1=0.15)
    # For analysing how does the amount of data (i.e. N) affects the estimation accuracy
    # and how does changing omega affect the estimation quality.
    # NOTE: N and omega are set based on the given simulation.name string!!!
    
    # select the true parameters etc. for generating the 'fake' data
    w <- regexpr('w',simulation.name)[1]
    r <- regexpr('rep',simulation.name)[1]
    true.params <- list()
    true.params$N <- as.integer(substr(simulation.name,19,w-2))
    true.params$true.neff.ind <- 22 # 22 is ~2000 in 100x100 grid
    true.params$true.mu.ind <- 22 # 22 is ~0.001 in 100x100 grid
    true.params$omega1 <- as.integer(substr(simulation.name,w+2,r-2)) / 100
    true.params$lambda <- 0.0001
    true.params$kt <- 5 # This is actually not estimated but assumed as fixed and known 
    true.params$hosp.times <- 100*c(5,10,15)
    
    opt3$use.prior.from.abc <- 0
    
    opt3$kt <- 5
    
    ## Gibbs sampling settings
    opt3$nr.chains <- 1
    opt3$nr.gibbs.samples <- 25000
    opt3$rand.gibbs.init <- TRUE # whether to initialize by drawing from the prior
    
    # some rough way of setting the proposal q
    if (true.params$N <= 50) {
      opt3$proposal.v.neff <- 750^2 # proposal variances in the Metropolis step
      opt3$proposal.v.mu <- 0.0003^2
    } else if (true.params$N <= 200) {
      opt3$proposal.v.neff <- 500^2 
      opt3$proposal.v.mu <- 0.0002^2
    } else {
      opt3$proposal.v.neff <- 200^2 
      opt3$proposal.v.mu <- 0.00015^2
    }
    
    opt3$which.data <- 'artificial1'
    opt3$which.arm <- NA
    opt3$abc.simulation.name <- 'abc_final_simul'
    opt3$cond.dens.simulation.name <- 'cond_dens_final_simul'
    
    opt3$alpha_lambda <- 2.5
    opt3$beta_lambda <- 1600
    opt3$kt <- 5
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  } else if (simulation.name == 'mixture_true_data_E_final' || simulation.name == 'mixture_true_data_E_final_pr') {
    
    ## Gibbs sampling settings
    opt3$nr.chains <- 4
    opt3$nr.gibbs.samples <- 25000
    opt3$rand.gibbs.init <- TRUE # whether to initialize by drawing from the prior
    
    opt3$proposal.v.neff <- 150^2 # proposal variances in the Metropolis step
    opt3$proposal.v.mu <- 0.0001^2
    
    opt3$use.prior.from.abc <- grepl('_pr', simulation.name)
    opt3$chosen.discrepancy.name <- 'l1'
    
    opt3$which.data <- 'true1'
    opt3$which.arm <- 'E'
    opt3$abc.simulation.name <- 'abc_final_simul'
    opt3$cond.dens.simulation.name <- 'cond_dens_final_simul'
    
    # these produce approximately d_mean=50=d_std
    #opt3$alpha_lambda <- 5
    #opt3$beta_lambda <- 7500
    #opt3$kt <- 5
    
    # these produce approximately d_mean=30, d_std=50
    opt3$alpha_lambda <- 2.5
    opt3$beta_lambda <- 1600
    opt3$kt <- 5
    
  } else {
    stop(paste('Incorrect simulation.name in task 3: ', simulation.name, sep = ''))
  }
  
  opt3$true.params <- true.params
  if (print.settings) { print(opt3) }
  return(opt3) 
}


####################################################################################

gamma.params <- function(m,v) {
  # Returns Gamma-distribution parameters alpha, beta from the given mean m and variance v.
  # Validity of inputs is not checked.
  alpha <- m^2/v
  beta <- m/v
  return(c(alpha, beta))
}

cg.params <- function(m,v,k=1) {
  # Returns Compound Gamma distribution parameters alpha, beta from the given mean m, variance v and k parameters.
  # If k == 1, then this is the Lomax distribution.
  # Validity of inputs is not fully checked, it is assumed that k >= 1.
  if (k*v < m^2) {
    stop('Invalid parameters for CG distribution.')
  }
  alpha <- (2*k*v + (k-1)*m^2)/(k*v - m^2)
  beta <- m*(v + m^2)/(k*v - m^2)
  return(c(alpha, beta))
}


####################################################################################

gen.neff.mu_grid <- function(opt=list(), grid.level='test') {
  # Generates 2d grid of n_eff,mu-values with different scarsity.
  # The grid and some related values are appended to list 'opt' that is required as input.
  
  if (grid.level == 'test') {
    # SIMPLE "1-POINT" GRID, ONLY FOR TESTING/DEBUGGING
    
    opt$n1 <- 1 # Grid length for n_eff (n.strains)
    opt$n2 <- 1 # Grid length for mu (mutation.rate)
    
    opt$neff.grid <- 10000
    opt$neff.min <- opt$neff.grid; opt$neff.max <- opt$neff.grid
    opt$mu.grid <- 0.0009
    opt$mu.min <- opt$mu.grid; opt$mu.max <- opt$mu.grid
    
  } else if(substr(grid.level,1,5) == 'final') {
    # GRID FOR ACTUAL COMPUTATIONS, MUST BE OF THE FORM 'final30x40' ETC.
    
    x <- regexpr('x',grid.level)[1]
    opt$n1 <- as.integer(substr(grid.level,6,x-1)) # Grid length for n_eff (n.strains)
    opt$n2 <- as.integer(substr(grid.level,x+1,nchar(grid.level))) # Grid length for mu (mutation.rate)
    if (is.na(opt$n1) || is.na(opt$n2)) {
      stop('Error with grid lengths.')
    }
    
    opt$neff.min <- 20
    opt$neff.max <- 10000
    opt$neff.grid <- get.int.grid(opt$neff.min, opt$neff.max, opt$n1, log.scale = FALSE)
    
    opt$mu.min <- 0.0002
    opt$mu.max <- 0.008
    opt$mu.grid <- seq(opt$mu.min, opt$mu.max, len=opt$n2)

  } else if(substr(grid.level,1,5) == 'paper') {
    # GRID FOR FINAL COMPUTATIONS FOR THE PAPER, MUST BE OF THE FORM 'paper30x40' ETC.
    
    x <- regexpr('x',grid.level)[1]
    opt$n1 <- as.integer(substr(grid.level,6,x-1)) # Grid length for n_eff (n.strains)
    opt$n2 <- as.integer(substr(grid.level,x+1,nchar(grid.level))) # Grid length for mu (mutation.rate)
    if (is.na(opt$n1) || is.na(opt$n2)) {
      stop('Error with grid lengths.')
    }
    
    opt$neff.min <- 20
    opt$neff.max <- 10000
    opt$neff.grid <- get.int.grid(opt$neff.min, opt$neff.max, opt$n1, log.scale = FALSE)
    
    opt$mu.min <- 0.00005
    opt$mu.max <- 0.005
    opt$mu.grid <- seq(opt$mu.min, opt$mu.max, len=opt$n2)
    
  } else {
    stop('Incorrect n_eff,mu-grid.')
  }
  opt$grid.name <- grid.level
  return(opt)
}





