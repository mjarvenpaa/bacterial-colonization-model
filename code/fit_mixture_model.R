fit.mixture.model <- function(root, simulation.name=NULL, fake.data.seed=123456) {
  # Uses Gibbs sampling and the ABC inference approximation (by earlier simulations) to estimate 
  # the parameters in the mixture model. 
  
  graphics.off()
  
  if (is.null(simulation.name)) {
    simulation.name <- 'mixture_test_artificial' # USE ARTIFICIAL DATA
  }
  
  ## General settings for fitting the mixture model
  plot.cond.dens.fig <- TRUE # whether plot the cond dens computed by simulations
  
  print.gibbs.init <- T # whether to print the initialisation point
  
  use.metropolis.step <- 1 # if true then Metropolis step and grid interpolation is used, otherwise grid sampling
  print.metropolis.info <- 0
  
  ###################################################################################################
  
  ## libraries, functions
  code.path <- paste(root, '/code', sep='')
  source(paste(code.path, '/get_settings.R', sep=''))
  source(paste(code.path, '/get_data.R', sep=''))
  source(paste(code.path, '/abc_methods_simul_model.R', sep=''))
  source(paste(code.path, '/cond_dens_help_functions.R', sep=''))
  source(paste(code.path, '/mixture_model.R', sep=''))
  source(paste(code.path, '/figure_functions.R', sep=''))
  source(paste(code.path, '/auxiliary_functions.R', sep=''))
  
  res.dir.name <- paste(root, '/simulation_outputs_mixture_model/', simulation.name, sep='')
  if (!file.exists(res.dir.name)) {
    dir.create(res.dir.name)
  }
  
  # Settings for mixture model and fitting
  opt3 <- get.mixture.fitting.settings(simulation.name, print.settings=TRUE)
  
  # Settings for getting cond dens
  opt2 <- get.cond.dens.settings(simulation.name=opt3$cond.dens.simulation.name, print.settings=FALSE)
  
  # Settings for prior of n_eff and mu which comes from earlier ABC inference
  opt <- get.abc.fit.settings(simulation.name=opt3$abc.simulation.name, print.settings=FALSE)
  
  # get some other settings
  use.prior.from.abc <- opt3$use.prior.from.abc
  chosen.discrepancy.name <- opt3$chosen.discrepancy.name
  proposal.v.neff <- opt3$proposal.v.neff # proposal variances in the Metropolis step
  proposal.v.mu <- opt3$proposal.v.mu
  print.output <- opt3$print.output
  
  ###################################################################################################
  
  ## Set up the grid of cond dens computation
  neff.grid.new <- opt2$neff.grid
  mu.grid.new <- opt2$mu.grid 
  
  ## grid for abc prior; this can be different grid than above now
  neff.grid.abc <- opt$neff.grid
  mu.grid.abc <- opt$mu.grid 
  
  ###################################################################################################
  
  ## Get and plot the CONDITIONAL DENSITIES computed previously in a n_eff,mu grid
  figure.filename <- paste(res.dir.name, '/p_di_densities', sep = '')
  figure.filename.nice.plot <- paste(res.dir.name, '/p_di_densities_nice_plot.pdf', sep = '')
  figure.filename.d.means <- paste(res.dir.name, '/p_di_densities_means.pdf', sep = '')
  p.sim.given.neff.mu <- get.p.sim.given.neff.mu.densities(root, opt2, neff.grid.new, mu.grid.new, figure.filename, 
                                                           figure.filename.nice.plot, figure.filename.d.means)

  ###################################################################################################
  
  # set seed for fake data generation (if true data is used, then this is redundant)
  set.seed(fake.data.seed)
  
  ## GET THE DATA
  clear.data <- get.mixture.model.data(root, opt3$which.data, opt3$which.arm, neff.grid.new, mu.grid.new, 
                                       p.sim = p.sim.given.neff.mu, plot.data = FALSE, print.info = FALSE, 
                                       true.params = opt3$true.params)
  d_is <- clear.data$d_is
  t_is <- clear.data$t_is
  td_metadata <- clear.data$td_metadata
  N <- length(d_is) # total number of pairs 
  
  ###################################################################################################
  
  ## Get the prior of n_eff,mu i.e. the posterior p(n_eff,mu|data_1:2)
  # normalisation is not necessary
  if (use.prior.from.abc) {
    abc.post.filename <- paste(root, '/simulation_outputs_abc/', opt3$abc.simulation.name, '/abc_results/', 
                               'ABC_samples_', chosen.discrepancy.name, '.RData', sep = '')
    load(file=abc.post.filename) 
    log.neff.mu.prior <- log(abc.samples$samples)
    # if no samples is observed < threshold, then the ABC posterior is 0. Fix issues caused by this
    # by setting it small but nonzero everywhere
    log.neff.mu.prior <- pmax(log.neff.mu.prior,log(0.01))
    
    # Since the prior density estimate can be a bit rough due to limited number of ABC samples, we use 2d kde 
    # to smooth it slightly 
    # NOTE: n_eff is discrete and samples are at grid points but this does not seem to matter here
    # NOTE2: this does not seem to matter much after all. In E data case the posterior does
    # look a bit unsmooth due to the grid approximation but it does not really affect the results
    smooth.prior <- FALSE
    if (smooth.prior) {
      library(MASS)
      totn <- sum(abc.samples$samples)
      x <- rep(NA,totn); y <- rep(NA,totn)
      c <- 1
      for (i in 1:length(neff.grid.abc)) {
        for (j in 1:length(mu.grid.abc)) {
          if (abc.samples$samples[i,j] > 0) {
            x[c:(c+abc.samples$samples[i,j]-1)] <- neff.grid.abc[i]
            y[c:(c+abc.samples$samples[i,j]-1)] <- mu.grid.abc[j]
            c <- c + abc.samples$samples[i,j]
          }
        }
      }
      kde.res <- kde2d(x = x, y = y, n = c(length(neff.grid.abc),length(mu.grid.abc)),
                       #h = c(100,0.0001), # bandwidths
                       lims = c(min(neff.grid.abc),max(neff.grid.abc),min(mu.grid.abc),max(mu.grid.abc)))
      log.neff.mu.prior <- log(kde.res$z)
      
      ## Plot the prior used for the inference
      if(1) {
        dev.new()
        # plot the samples
        contour(abc.samples$samples, levels = c(1,2,5,10,50,100,1000))
        grid()
        dev.new()
        # plot the estimated density by 2d kde
        contour(exp(log.neff.mu.prior), levels = c(0.001,0.01,0.1,0.2,0.3,0.4,0.5,1))
        grid()
      }
    }
    
  } else {
    # just use uniform prior
    log.neff.mu.prior <- matrix(0, length(neff.grid.abc), length(mu.grid.abc))
  }

  
  ###################################################################################################
  ###################################################################################################
  
  ## Get model related settings
  kt <- opt3$kt
  
  ## Get Gibbs sampling settings
  nr.chains <- opt3$nr.chains
  nr.gibbs.samples <- opt3$nr.gibbs.samples
  rand.gibbs.init <- opt3$rand.gibbs.init # if true, initial points are drawn from the prior
  
  ## Get (Hyper)priors 
  alpha_lambda <- opt3$alpha_lambda
  beta_lambda <- opt3$beta_lambda
  gamma_omega <- opt3$gamma_omega
  
  ## Other (debugging) settings etc.
  print.freq <- floor(ifelse(nr.gibbs.samples<=1000,nr.gibbs.samples/10,nr.gibbs.samples/20))
  neff.mu.debug.plot <- TRUE
  neff.mu.debug.plot.freq <- 200
  plot.indiv.chains <- TRUE
  z.plot.inds <- get.z.plot.inds(d_is)
  
  #################################################################################################
  
  ## PRECOMPUTE p_S DENSITIES
  log.ps.dens <- array(-Inf, c(length(neff.grid.new), length(mu.grid.new), N))
  for (i in 1:length(neff.grid.new)) {
    for (j in 1:length(mu.grid.new)) {
      if (any(is.na(p.sim.given.neff.mu[[i]][[j]]))) {
        # All simulations missing due to prior bounds
        next
      }
      log.ps.dens[i,j,] <- log.p_S.density(d_is, t_is, mu.grid.new[j], mu.grid.new, p.sim.given.neff.mu[[i]][[j]])
    }
  }
  
  #################################################################################################
  
  ## DIFFERENT CHAINS
  for (cur.chain in 1:nr.chains) { 
    
    ## Set different seed for each chain
    set.seed(123456 + cur.chain) 
    
    ## Initialize variables
    zs <- array(NA,c(N,2,nr.gibbs.samples+1)) # N x 2 x nr.gibbs.samples, each row is (z_i1,z_i2)
    omegas <- matrix(NA,2,nr.gibbs.samples+1) # 2 x nr.gibbs.samples
    lambdas <- rep(NA,nr.gibbs.samples+1) # 1 x nr.gibbs.samples
    neffs <- rep(NA,nr.gibbs.samples+1) # 1 x nr.gibbs.samples
    mus <- rep(NA,nr.gibbs.samples+1) # 1 x nr.gibbs.samples
    t_0s <- matrix(NA,N,nr.gibbs.samples+1) # N x nr.gibbs.samples
    eettas <- matrix(NA,N,nr.gibbs.samples+1) # N x nr.gibbs.samples
    log.unnorm.post <- rep(NA,nr.gibbs.samples+1) # unnormalised posterior values in the chain
    
    ## Set initial values
    if (rand.gibbs.init) {
      ## Random initialization from the prior
      prior.gen <- rand.from.prior(N, alpha_lambda, beta_lambda, gamma_omega, neff.grid.abc, 
                                   mu.grid.abc, kt, log.neff.mu.prior)
      zs[,,1] <- prior.gen$z
      omegas[,1] <- prior.gen$omega
      lambdas[1] <- prior.gen$lambda
      neffs[1] <- prior.gen$neff
      mus[1] <- prior.gen$mu
      t_0s[,1] <- prior.gen$t_0
      eettas[,1] <- prior.gen$mu * prior.gen$t_0
      
      # print selected initial location
      if (print.gibbs.init) {
        print('Initialisation from the prior:')
        prior.gen$z <- t(prior.gen$z) # for a bit nicer printing
        print(prior.gen)
      }
    } else {
      ## Fixed initial point
      zs[,,1] <- cbind(rep(FALSE,N),rep(TRUE,N))
      omegas[,1] <- c(0.75,0.25)
      lambdas[1] <- 0.001
      neffs[1] <- neff.grid.new[1]
      mus[1] <- mu.grid.new[1]
      t_0s[,1] <- rep(1000,N)
      eettas[,1] <- mus[1] * t_0s[,1]
      
      if (print.gibbs.init) {
        print('Using fixed initial point for Gibbs sampling.')
      }
    }
    
    if (print.output && nr.chains > 1) {
      print(paste('chain=',cur.chain,'/',nr.chains,sep=''))
    }
    neff.mu.acc.moves <- 0
    
    #################################################################################################
    
    ## Run (systematic scan) Gibbs sampling 
    for (i in 2:(nr.gibbs.samples+1)) {
      
      if (print.output && (i-1 <= 5 || !(i-1)%%print.freq)) {
        print(paste('Iteration: ',i-1,'/',nr.gibbs.samples,sep = ''))
        if (i >= 100 && use.metropolis.step) {
          print(paste('Acceptance prob. = ', neff.mu.acc.moves/(i-1), sep=''))
        }
      }
      
      ## 1/5: simulate neff,mu jointly (using the 2d grid values)
      ###########################################################
      
      sum.zt <- sum(zs[,2,i-1] * t_is)
      sum.eettas <- sum(eettas[,i-1])
      
      #===================================================================================================
      if (use.metropolis.step) {
        ## Use Metropolis step and interpolation
        cur.neff.mu <- c(neffs[i-1], mus[i-1])
        
        # propose new point
        proposed.neff.mu <- gen.neff.mu.proposal(cur.neff.mu, proposal.v.neff, proposal.v.mu)
        
        # check if we are in prior support
        inside.bounds <- prior.neff.mu.value(proposed.neff.mu[1], proposed.neff.mu[2], 
                                             neff.grid.new, mu.grid.new, check.grid=FALSE)
        if (inside.bounds) {
          
          # new proposed point is inside the prior bounds
          # compute the conditional probability both at the old point and the proposed point 
          log.neff.mu.cond.dens.interp <- rep(NA,2)
          for (j in 1:2) {
            # compute first old point, then the new point
            neff.mu.j <- (j==1)*cur.neff.mu + (j==2)*proposed.neff.mu
            pt.inds <- get.2x2.grid(neff.grid.new, mu.grid.new, neff.mu.j)
            n1j <- pt.inds[1]
            m1j <- pt.inds[2]
            
            zs.log.ps.dens <- sweep(log.ps.dens[n1j:(n1j+1),m1j:(m1j+1),], 3, zs[,1,i-1], '*') 
            sum.zs.log.ps.dens <- apply(zs.log.ps.dens, 1:2, sum) 
            
            # add log prior
            #sum.zs.log.ps.dens <- sum.zs.log.ps.dens + log.neff.mu.prior[n1j:(n1j+1),m1j:(m1j+1)] 
            
            # interpolate
            #sum.zs.log.ps.dens[!is.finite(sum.zs.log.ps.dens)] <- -1e-20
            sum.zs.log.ps.dens.interp <- (interp.bilinear.vectorized((sum.zs.log.ps.dens), neff.grid.new[n1j:(n1j+1)], 
                                                                    mu.grid.new[m1j:(m1j+1)], 1, 1, neff.mu.j))
            
            # compute prior density value by interpolation
            pt.inds.abc <- get.2x2.grid(neff.grid.abc, mu.grid.abc, neff.mu.j)
            n1j.abc <- pt.inds.abc[1]
            m1j.abc <- pt.inds.abc[2]
            prior.interp <- interp.bilinear.vectorized(log.neff.mu.prior[n1j.abc:(n1j.abc+1),m1j.abc:(m1j.abc+1)], 
                                                       neff.grid.abc[n1j.abc:(n1j.abc+1)], 
                                                       mu.grid.abc[m1j.abc:(m1j.abc+1)], 1, 1, neff.mu.j)
            
            # compute the rest of the terms
            nmtij <- 2*eettas[,i-1] + neff.mu.j[2]*t_is
            mt1j <- sum(d_is * zs[,2,i-1]*log(nmtij)) - neff.mu.j[2] * sum.zt
            #mt2j <- -lambdas[i-1] / neff.mu.j[2] * sum.eettas # if conditioned on lambda
            mt2j <- -(N*kt + alpha_lambda) * log(sum.eettas/neff.mu.j[2] + beta_lambda) # if lambda is collapsed 
            mt3j <- -N*kt*log(neff.mu.j[2])
            sum.term <- mt1j + mt2j + mt3j
            log.neff.mu.cond.dens.interp.j <- sum.zs.log.ps.dens.interp + sum.term + prior.interp
            
            #log.neff.mu.cond.dens.interp.j[is.nan(log.neff.mu.cond.dens.interp.j)] <- -Inf
            log.neff.mu.cond.dens.interp[j] <- log.neff.mu.cond.dens.interp.j
          }
          # compute the ratio
          log.roo <- (log.neff.mu.cond.dens.interp[2] - log.neff.mu.cond.dens.interp[1])
        
        } else {
          
          # outside bounds, keep the old point
          log.roo <- -Inf
          accept <- FALSE
          if (print.metropolis.info) {
            print('Proposed point was outside prior bounds.')
          }
        }
        
        if (is.null(log.roo)) {
          warning('Acceptance check results NULL.', immediate. = TRUE) # should not happen
          browser()
          
        } else if (is.na(log.roo)) {
          warning('Acceptance check results NaN/NA.', immediate. = TRUE)
          if (is.finite(log.neff.mu.cond.dens.interp[2])) {
            # old (initial) point was bad but the new one appears to be ok, so accept it and try to go on
            accept <- TRUE
            #browser()
          } else if (is.finite(log.neff.mu.cond.dens.interp[1])) {
            # old (initial) point is ok but the new one is bad, so reject it and try to go on
            accept <- FALSE
            #browser()
          } else {
            # both points are bad >:(
            # ad hoc resolution: select new point that should be ok and try to go on
            proposed.neff.mu <- c(neff.grid.new[1],mu.grid.new[1])
            accept <- TRUE
            #browser()
          }
          
        } else if (is.infinite(log.roo) && log.roo > 0) {
          #warning('Acceptance check results Inf.', immediate. = TRUE) 
          # old point was -Inf and new point is any finite value so we can accept it 
          # (with probability one) and try to go on
          accept <- TRUE
        } else {
          # everything was ok and a proper value for the acceptance criterion was computed
          accept <- (log.roo >= log(runif(n = 1)))
        }
        
        # check acceptance
        if (accept) {
          neffs[i] <- proposed.neff.mu[1]
          mus[i] <- proposed.neff.mu[2]
          neff.mu.acc.moves <- neff.mu.acc.moves + 1
          if (print.metropolis.info) { print('Accepted') }

        } else {
          neffs[i] <- cur.neff.mu[1]
          mus[i] <- cur.neff.mu[2]
          if (print.metropolis.info) { print('Rejected') }

        }
        #===================================================================================================
        
      } else {
        ## THIS IS AN ALTERNATIVE (BUT MORE SLOW) SAMPLING METHOD!
        # NOTE: WE ASSUME THAT THE GRID FOR COND.DENS AND ABC PRIOR COMPUTATIONS ARE THE SAME HERE
        
        ## Use the conditional density computed and sampled in the grid
        zs.log.ps.dens <- sweep(log.ps.dens, 3, zs[,1,i-1], '*')
        sum.zs.log.ps.dens <- colSums(aperm(zs.log.ps.dens, c(3,1,2)))
        
        muti <- outer(t_is, mu.grid.new)
        mt1j <- colSums((d_is*zs[,2,i-1])*(log(2*eettas[,i-1] + muti)))
        
        mt2j <- -mu.grid.new * sum.zt
        #mt3j <- -lambdas[i-1] / mu.grid.new * sum.eettas # if conditioned on lambda
        mt3j <- -(N*kt + alpha_lambda) * log(sum.eettas/mu.grid.new + beta_lambda) # if lambda is collapsed 
        mt4j <- -N*kt*log(mu.grid.new)
        sum.term <- mt1j + mt2j + mt3j + mt4j
        log.neff.mu.cond.dens <- t(t(sum.zs.log.ps.dens) + sum.term)

        if (!all(is.finite(log.neff.mu.cond.dens))) {
          # shouldn't happen anymore...
          #stop('NA/Inf problem with log.neff.mu.cond.dens.')

          # NaN can happen because prior bounds... set them to -Inf
          log.neff.mu.cond.dens[is.nan(log.neff.mu.cond.dens)] <- -Inf
        }
        # add log prior
        log.neff.mu.cond.dens <- log.neff.mu.cond.dens + log.neff.mu.prior

        neff.mu.inds <- sample.cat.density.matrix(log.neff.mu.cond.dens)
        r.neff.ind <- neff.mu.inds[1]
        r.mu.ind <- neff.mu.inds[2]
        neffs[i] <- neff.grid.new[r.neff.ind] 
        mus[i] <- mu.grid.new[r.mu.ind]


        ## Debugging: plot the conditional density for n_eff,mu
        if (neff.mu.debug.plot && (i==2 || (i-1)%%neff.mu.debug.plot.freq == 0)) {

          # compute (non-log) bin probabilities by normalising (for debug plotting only)
          # NOTE: the n_eff grid may not actually be exactly even due to rounding to integers
          neff.mu.cond.dens <- exp(log.neff.mu.cond.dens - log.sum.exp(log.neff.mu.cond.dens))
          neff.mu.cond.dens <- neff.mu.cond.dens / ((mu.grid.new[2] - mu.grid.new[1]) * sum(neff.mu.cond.dens))

          if (i==2) {
            require(grDevices) # for colours
            dev.new()
          }
          if (all(is.finite(neff.mu.cond.dens)) && any(neff.mu.cond.dens > 0)) {
            # plot density if numerics have not failed (due to underflow etc)
            titl <- paste('neff, mu cond. density, chain=',cur.chain,', iteration=',i-1,sep='')
            filled.contour(x=neff.grid.new, y=mu.grid.new, z=neff.mu.cond.dens, nlevels = 50, color = terrain.colors,
                           plot.title = title(titl,xlab='n_eff',ylab='mu'))
          } else {
            # plot log density otherwise
            titl <- paste('neff, mu log cond., chain=',cur.chain,', iteration=',i-1,sep='')
            filled.contour(x=neff.grid.new, y=mu.grid.new, z=log.neff.mu.cond.dens, nlevels = 50, color = terrain.colors,
                           plot.title = title(titl,xlab='n_eff',ylab='mu'))
          }
        }
      }
      
      
      ## 2/5: simulate lambda (analytic)
      ##################################
      if (!use.metropolis.step || (use.metropolis.step && accept)) {
        # NOTE: IF MU IS NOT INTEGRATED OUT IN METROPOLIS STEP, THEN NEW LAMBDA SHOULD ALWAYS BE GENERATED HERE 
        # under the grid computation approximation, we sample exactly and we always accept 
        lambdas[i] <- rgamma(n = 1, shape = N*kt + alpha_lambda, rate = sum.eettas/mus[i] + beta_lambda)
      } else {
        lambdas[i] <- lambdas[i-1]
      }
      
      
      ## 3/5: simulate eetta_i / t_0i (analytic)
      ##########################################
      ind.gam.cases <- which(zs[,1,i-1])
      ind.gam2.cases <- which(zs[,2,i-1] & d_is == 0)
      ind.gm.cases <- which(zs[,2,i-1] & d_is > 0)
      
      # Gamma/Exp cases, easy sampling
      eettas[ind.gam.cases,i] <- rgamma(n = length(ind.gam.cases), shape = kt, rate = lambdas[i]/mus[i])
      eettas[ind.gam2.cases,i] <- rgamma(n = length(ind.gam2.cases), shape = kt, rate = (2 + lambdas[i]/mus[i]))
      
      # gamma mixture case
      eettas[ind.gm.cases,i] <- gen.eettas(di = d_is[ind.gm.cases], ti = t_is[ind.gm.cases], mu = mus[i], 
                                           lambda = lambdas[i], kt = kt)
      
      # keep track of t_0s
      t_0s[,i] <- eettas[,i] / mus[i]
      
      
      ## 4/5: simulate z (analytic, but depends on p_sim)
      ###################################################
      if (use.metropolis.step) {
        # Metropolis case, here we need to use interpolation, again
        if (!exists("log.ps.dens.interp") || neffs[i] != neffs[i-1] || mus[i] != mus[i-1]) {
          # compute only if the point has changed or the first time here
          cur.neff.mu <- c(neffs[i],mus[i])
          cur.bd.inds <- get.2x2.grid(neff.grid.new, mu.grid.new, cur.neff.mu)
          log.ps.dens.interp <- interp.bilinear.vectorized(log.ps.dens, neff.grid.new, mu.grid.new, 
                                                           cur.bd.inds[1], cur.bd.inds[2], cur.neff.mu)
        }
        log.pzi1 <- log(omegas[1,i-1]) + log.ps.dens.interp
      } else {
        # grid computation case
        log.pzi1 <- log(omegas[1,i-1]) + log.ps.dens[r.neff.ind, r.mu.ind,]
      }
      
      nmti <- 2*eettas[,i] + mus[i]*t_is
      log.pzi2 <- log(omegas[2,i-1]) + (d_is * log(nmti) - nmti - lfactorial(d_is))
      
      r <- sample.cat.density.vectorized(cbind(log.pzi1,log.pzi2))
      zs[,1,i] <- !as.logical(r-1)
      zs[,2,i] <- !zs[,1,i]
      
      if (any(!is.finite(zs[,1,i]))) {
        warning('Something wrong with z!', immediate. = TRUE)
        browser()
      }
      
      
      ## 5/5: simulate omega (analytic)
      #################################
      # Using Beta simulation for sampling from 2d Dirichlet
      n1 <- sum(zs[,1,i])
      n2 <- sum(zs[,2,i])
      omegas[1,i] <- rbeta(n = 1, shape1 = n1 + gamma_omega[1], shape2 = n2 + gamma_omega[2])
      omegas[2,i] <- 1 - omegas[1,i]
      
      
      ## Compute unnormalised log posterior (for plotting/convergence inspection) 
      ## (analytic, but depends on p_sim and prior of (n_eff,mu))
      ###########################################################################
      if (use.metropolis.step) {
        zs.log.ps.dens <- zs[,1,i] * log.ps.dens.interp
        if (!exists("prior.term") || neffs[i] != neffs[i-1] || mus[i] != mus[i-1]) {
          cur.neff.mu <- c(neffs[i],mus[i])
          cur.bd.inds.abc <- get.2x2.grid(neff.grid.abc, mu.grid.abc, cur.neff.mu)
          prior.term <- interp.bilinear.vectorized(log.neff.mu.prior, neff.grid.abc, mu.grid.abc,
                                                   cur.bd.inds.abc[1], cur.bd.inds.abc[2], cur.neff.mu)
        }
      } else {
        # grid computation case
        zs.log.ps.dens <- zs[,1,i] * log.ps.dens[r.neff.ind, r.mu.ind,]
        prior.term <- log.neff.mu.prior[r.neff.ind, r.mu.ind]
      }
      
      sum.eettas <- sum(eettas[,i])
      term11 <- sum(zs.log.ps.dens)
      term12 <- sum(zs[,2,i] * (d_is*log(nmti) - nmti - lfactorial(d_is)))
       
      term2 <- (n1 + gamma_omega[1] - 1)*log(omegas[1,i]) + (n2 + gamma_omega[2] - 1)*log(omegas[2,i])
      term3 <- (N*kt + alpha_lambda - 1)*log(lambdas[i]) - (sum.eettas/mus[i] + beta_lambda)*lambdas[i]
      termn <- (kt - 1)*sum(log(eettas[,i])) - N*kt*log(mus[i])
      log.unnorm.post[i] <- term11 + term12 + term2 + term3 + termn + prior.term
    }
    
    if (use.metropolis.step) {
      print(paste('Final acceptance prob. = ', neff.mu.acc.moves/nr.gibbs.samples, sep=''))
    }
    
    ## Save computed results
    sample.filename <- paste(res.dir.name, '/mixture_fitting_test', cur.chain, '.RData', sep='')
    if (cur.chain == 1) {
      # save clear data etc. for possible later inspections only once
      # do not save eettas to save disk space and because they can be recovered from t_0 and mu
      save(clear.data, d_is, t_is, zs, omegas, lambdas, neffs, mus, t_0s, log.unnorm.post, cur.chain, 
           file=sample.filename)
    } else {
      save(zs, omegas, lambdas, neffs, mus, t_0s, log.unnorm.post, cur.chain, file=sample.filename)
    }
    
    ## Plot and print results in current chain 
    if (plot.indiv.chains) {
      sample.filename <- paste(res.dir.name, '/fig_samples', cur.chain, '.pdf', sep='')
      simple.sample.plot(zs, omegas, lambdas, neffs, mus, t_0s, log.unnorm.post, t_is, d_is, 
                         neff.grid.new, mu.grid.new, z.plot.inds, sample.filename) 
      
      sample.filename2 <- paste(res.dir.name, '/fig_neff_vs_mu', cur.chain, '.pdf', sep='')
      plot.neff.mu.joint.post(neffs, mus, neff.grid.new, mu.grid.new, fig.name=sample.filename2, 
                              plot.density.kde = use.metropolis.step, true.params=clear.data$true.params)
    }
  }
}



  