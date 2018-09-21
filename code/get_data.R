# Functions for reading the data and returning it in a nice format. 


get.mixture.model.data <- function(root, which.data='true1', which.arm='E', neff.grid=NULL, mu.grid=NULL, p.sim = NULL, plot.data=TRUE, 
                                   true.params = NULL, print.info = FALSE) {
  # Wrapper for loading the data needed for fitting the mixture model
  
  ## Select the data that is returned and used for model fitting
  if (which.data == 'artificial1') {
    warning('Artificial data returned for testing.', immediate. = TRUE)
    return(get.artificial.mixture.data(neff.grid, mu.grid, p.sim, true.params, print.info))
    
  } else if (which.data == 'true1') {
    #return(get.mixture.model.data.true(root, arm = which.arm, plot.data = plot.data, print.info=print.info, which.data = which.data))
    return(get.preprocessed.mrsa.data(root, which.arm))
  } else {
    stop('Incorrect data set requested.')
  }
}


get.artificial.mixture.data <- function(neff.grid, mu.grid, p.sim, true.params, print.info) {
  # Generates randomly data from the model to be used for testing/demonstrating that 
  # the estimation algorithm itself works.
  # Settings on how to generate the data are given as input in 'true.params'
  
  true.params$neff <- neff.grid[true.params$true.neff.ind] 
  true.params$mu <- mu.grid[true.params$true.mu.ind] 
  
  # generate data
  true.params$t0i <- rgamma(n = true.params$N, shape = true.params$kt, rate = true.params$lambda) 
  probs <- c(true.params$omega1, 1 - true.params$omega1)
  true.params$z1 <- sample(x = c(TRUE,FALSE), size = true.params$N, replace = TRUE, prob = probs)
  nm1 <- length(which(true.params$z1))
  nm2 <- length(which(!true.params$z1))
  
  ti <- sample(x = true.params$hosp.times, size = true.params$N, replace = TRUE)
  di <- rep(NA,true.params$N)
  
  # from model1:
  # generate from the simulated cond dens
  probs <- p.sim[[true.params$true.neff.ind]][[true.params$true.mu.ind]]
  di1 <- sample(x = 0:(length(probs)-1), size = nm1, replace = TRUE, prob = probs)
  di2 <- rpois(n = nm1, lambda = true.params$mu * ti[true.params$z1])
  di[true.params$z1] <- di1 + di2
  
  # from model2:
  mut <- true.params$mu * (2*true.params$t0i[!true.params$z1] + ti[!true.params$z1])
  di[!true.params$z1] <- rpois(n = nm2, lambda = mut)
  
  # order data
  ind <- order(di, decreasing = FALSE)
  di <- di[ind]
  ti <- ti[ind]
  true.params$t0i <- true.params$t0i[ind]
  true.params$z1 <- true.params$z1[ind]
  
  # print the generated data
  if (print.info) {
    print('--------------------')
    print('Artificial data:')
    print(paste('n_eff: ', true.params$neff, ', mu: ', true.params$mu, sep=''))
    print('t_i:'); print(ti)
    print('t_0i:'); print(true.params$t0i)
    print('z_1:'); print(true.params$z1)
    print('d_i:'); print(di)
    print('--------------------')
  }
  return(list(t_is = ti, d_is = di, td_metadata = NA, 
              t_is_excl = NA, d_is_excl = NA, td_metadata_excl = NA, true.params = true.params))
}


get.preprocessed.mrsa.data <- function(root, which.arm) {
  # Returns the t_i,d_i-pairs of the MRSA data that has been preprocessed and saved to a file. 
  
  # Loads the preprocessed datafile and returns its contents.
  # Note: 'metadata', which is necessary for linking the distances to individual patients, could not be 
  # included to preprocessed data and is thus not returned.
  
  data.filename <- 'mrsa_clear_testdata_arm'
  filename <- paste(root,'/data_CLEAR_minimal/',data.filename,which.arm,'.RData',sep = '')
  load(file = filename)
  
  return(list(t_is = t_is, d_is = d_is, td_metadata = NA, 
              t_is_excl = t_is_excl, d_is_excl = d_is_excl, td_metadata_excl = NA))
}


#####################################################################################################


get.sim.model.data <- function(root) {
  # Note: This function returns artificial data that is similar to real, observed data D_0 BUT that is 
  # generated from the model itself (using the function below). 
  
  data.filename <- 'distance_test_data'
  filename <- paste(root,'/data_external/',data.filename,'.RData',sep = '')
  load(file = filename)
  return(data)
}


get.artificial.abc.data <- function(root) {
  # Generate artificial test data that is similar to data D_0 used in the analysis.
  # Note: This function generates also data of all patients "1-8" although in the true data analysis
  # we actually neglect one of the patients ("patient #1219")
  
  set.seed(123456)
  
  code.path <- paste(root, '/code', sep='')
  source(paste(code.path, '/simulation_functions.R', sep=''))
  source(paste(code.path, '/data_summary_functions.R', sep=''))
  source(paste(code.path, '/abc_methods_simul_model.R', sep=''))
  
  # Set the parameter to use:
  # use the final estimates based on the full analysis with the true data:
  neff.true <- 1700
  mu.true <- 0.00076

  
  # these are as in true data D_0
  data <- list(data.first=vector('list',8),data.last=vector('list',8),data.AM=vector('list',13))
  p18.pairs <- c(66,rep(28,7))
  pAM.pairs <- c(45,45,55,66,45,36,45,45,36,66,45,45,36)
  p18.nr.genomes <- c(12,rep(8,7))
  pAM.nr.genomes <- c(10,10,11,12,10,9,10,10,9,12,10,10,9)
  
  genome.length <- 3000000
  
  # patients 1-8
  for (i in 1:8) {
    save.interval = c(1000,6000)
    pops.i <- simulate.evolution(mutation.rate=mu.true, n.strains=neff.true, genome.length=genome.length, 
                                    n.generations.to.simulate=save.interval[length(save.interval)], save.interval=save.interval, 
                                    save.path=NULL, save.extension=NULL, print.output=TRUE, save.output=FALSE)
    
    pops.i1 <- resample.pop(pops=pops.i, n.strains.new=p18.nr.genomes[i])
    data$data.first[[i]] <- distance.vector(compute.distance.distribution(pops.i1[[1]]),p18.pairs[i])
    names(data$data.first)[i] <- paste('patient',i,sep='')
    data$data.last[[i]] <- distance.vector(compute.distance.distribution(pops.i1[[2]]),p18.pairs[i])
    names(data$data.last)[i] <- paste('patient',i,sep='')
  }
  
  # patients A-M
  for (i in 1:13) {
    save.interval <- c(1,ceiling(runif(n=1,min=500,max=9000))) # generate randomly as assumed in the real data case
    pops.i <- simulate.evolution(mutation.rate=mu.true, n.strains=neff.true, genome.length=genome.length, 
                                 n.generations.to.simulate=save.interval[length(save.interval)], save.interval=save.interval, 
                                 save.path=NULL, save.extension=NULL, print.output=TRUE, save.output=FALSE)
    
    pops.i1 <- resample.pop(pops=pops.i, n.strains.new=pAM.nr.genomes[i])
    data$data.AM[[i]] <- distance.vector(compute.distance.distribution(pops.i1[[2]]),pAM.pairs[i])
    names(data$data.AM)[i] <- paste('patient',LETTERS[i],sep='')
    #names(data$data.AM)[i] <- paste('patient',LETTERS[i],' ',save.interval[2],sep='')
  }
  
  # print data
  #print(data)
  
  # save data to file
  data.filename <- 'distance_test_data'
  filename <- paste(root,'/data_external/',data.filename,'.RData',sep = '')
  save(data, file = filename)
  invisible(data)
}


distance.vector <- function(dist,supposed.len) {
  # Extracts the full distance vector from the distance datastructure.
  n <- dim(dist$snp.set.dist)[1]
  v <- NULL
  for (i in 1:n) {
    for (j in i:n) {
      v <- c(v,rep(dist$snp.set.dist[i,j], dist$dist.elem.count[i,j]))
    }
  }
  if (length(v) != supposed.len) {
    # check just in case that the distance vector has correct length
    stop('Length of distance vector is wrong.')
  }
  return(sort(v))
}




