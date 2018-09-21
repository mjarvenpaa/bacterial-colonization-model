# Functions for ABC-inference (and related visualisations) for the simulation model

abc.posterior.table.from.triton.simuls <- function(root) {
  # Loads the simulated discrepancies and computes the ABC posterior density in the grid.
  # Some visualisations are also made. 
  #
  # NOTE: THIS IS FOR LOCAL RUN. 
  # COMPUTING THE ABC POSTERIOR AND RELATED VISUALISATIONS IS CHEAP AFTER THE COMPUTATIONALLY COSTLY 
  # SIMULATIONS HAVE BEEN COMPLETED IN CLUSTER SO THIS CAN BE RUN LOCALLY.
  
  simulation.name <- 'abc_final_simul' 
  
  methods.to.run <- c('visualisation_discrepancy','ABC','visualisation_ABC')

  # which discrepancies to use for ABC posterior computations/plottings
  chosen.discrepancy.names <- list('eucl_noAM','eucl','l1_noAM','l1')
  
  ## libraries, functions
  code.path <- paste(root, '/code', sep='')
  source(paste(code.path, '/get_settings.R', sep=''))
  source(paste(code.path, '/get_data.R', sep=''))
  source(paste(code.path, '/auxiliary_functions.R', sep=''))
  
  res.dir.name <- paste(root, '/simulation_outputs_abc/', simulation.name, sep='')
  
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  ## Get the settings used in discrepancy simulations
  opt <- get.abc.fit.settings(simulation.name=simulation.name)
  
  parallel.runs <- opt$parallel.runs
  
  n.rep <- opt$n.rep
  print.output <- opt$print.output
  
  neff.grid <- opt$neff.grid
  mu.grid <- opt$mu.grid
  
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  # GO THROUGH THE DISCREPANCIES TO BE COMPUTED
  for (cur.discr in chosen.discrepancy.names) {
    ## Load the (parameter, discrepancy)-pairs from the files and add them to a single datastructure
    discr.data <- array(NA, c(length(neff.grid), length(mu.grid), n.rep))
    
    for (id.i in 1:parallel.runs) {
      if (print.output) {
        print(paste(id.i,'/',parallel.runs,sep = ''))
      }
      
      # load the discrepancy data from files
      discrepancy.filename.id <- paste(res.dir.name, '/discrepancies/', cur.discr, '_id', id.i, '.RData', sep='')
      load(discrepancy.filename.id) # discrepancy.data
      
      reps <- get.current.ids(id.i, parallel.runs, n.rep)
      discr.data[,,reps] <- discrepancy.data
    }
    
    # discr.data contains now the discrepancy values [neff,mu,rep] -> discrepancy
    if (print.output) {
      print('Discrepancy data loaded.')
      print(discr.data)
    }
    
    
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    ## Plot the discrepancy values as 2d plot
    if (any(methods.to.run == 'visualisation_discrepancy')) {
      
      figure.filename <- paste(res.dir.name, '/abc_results/', cur.discr, '_fig', '.pdf', sep = '')
      plot.discrepancy.table(discr.data, neff.grid, mu.grid, cur.discr, figure.filename, print.output)
    }
    
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    ## Compute the ABC posterior in 2d grid for mixture model fitting
    if (any(methods.to.run == 'ABC')) {
      
      # compute the ABC density based on the simulations
      abc.samples <- return.abc.post.samples(discr.data, opt) 
      
      abc.post.filename <- paste(res.dir.name, '/abc_results/', 'ABC_samples_', cur.discr, 
                                 #'_eps', signif(abc.samples$epsilon,6), 
                                 '.RData', sep = '')
      save(abc.samples, file=abc.post.filename)
    }
    
    ## Plot the 2d ABC posterior computed above
    if (any(methods.to.run == 'visualisation_ABC')) {
      
      titl <- paste('ABC posterior, "', cur.discr, '" discrepancy', sep = '')
      figure.filename <- paste(res.dir.name, '/abc_results/', cur.discr, 
                               #'_eps', signif(abc.samples$epsilon,6), 
                               '_fig_ABC', '.pdf', sep = '')
      plot.abc.table(abc.samples$samples, neff.grid, mu.grid, titl, figure.filename, print.output)
    }
  }
}

###############################################################################################

plot.abc.post.nice <- function(root, simulation.name) {
  # Plot some nice figure with multiple ABC results next to each other for the paper. 
  # One needs to run the above script first. 
  
  # which discrepancies/ABC to plot nicely
  chosen.discrepancy.names <- list('l1','l1_noAM')
  titles <- list('A','B')
  
  lvls <- seq(0,0.02,len=40) # adjust this according to the posteriors to be plotted
  
  ## libraries, functions
  code.path <- paste(root, '/code', sep='')
  source(paste(code.path, '/get_settings.R', sep=''))
  source(paste(code.path, '/get_data.R', sep=''))
  source(paste(code.path, '/auxiliary_functions.R', sep=''))
  res.dir.name <- paste(root, '/simulation_outputs_abc/', simulation.name, sep='')
  
  opt <- get.abc.fit.settings(simulation.name=simulation.name)
  neff.grid <- opt$neff.grid
  mu.grid <- opt$mu.grid
  
  figure.filename <- paste(res.dir.name, '/abc_results/', 'abc_comparison', '.pdf', sep = '')
  
  n.plots <- length(chosen.discrepancy.names)
  pdf(file=figure.filename, width = 4*n.plots, height = 4, pointsize = 12, bg='white')
  plot.new()
  #par(mfrow=c(1,n.plots))
  
  # For plotting only part of the parameter space:
  smaller.grid <- 1
  if (smaller.grid) {
    n <- c(20,10000) # plot limits for n_eff
    m <- c(0.00005,0.0025) # plot limits for mu
    ln <- max(which(neff.grid<=n[1])); un <- min(which(neff.grid>=n[2]))
    lm <- max(which(mu.grid<=m[1])); um <- min(which(mu.grid>=m[2]))
    neff.grid <- neff.grid[ln:un]
    mu.grid <- mu.grid[lm:um]
  } else {
    ln <- 1; un <- length(neff.grid)
    lm <- 1; um <- length(mu.grid)
  }
  
  i <- 1
  for (cur.discr in chosen.discrepancy.names) {
    # load data
    abc.post.filename <- paste(res.dir.name, '/abc_results/', 'ABC_samples_', cur.discr, 
                               '.RData', sep = '')
    load(file=abc.post.filename)
    samples <- abc.samples$samples[ln:un,lm:um]
    
    if (i == 1) {
      par(new = "TRUE", plt = c(0.075,0.45,0.15,0.865), mgp = c(1.7,0.55,0))
    } else if (i == 2) {
      par(new = "TRUE", plt = c(0.575,0.95,0.15,0.865), mgp = c(1.7,0.55,0))
    } else {
      # todo
      break
    }
    plot.abc.table(samples, neff.grid, mu.grid, titl = titles[[i]], figure.filename = NULL, 
                   plot.filled.contour = TRUE, lvls = lvls)
    i <- i + 1
  }
  dev.off()
  invisible(0)
}


#########################################################################################
#########################################################################################

plot.discrepancy.table <- function(discr.data, neff.grid, mu.grid, discrepancy.type, figure.filename, print.output=FALSE) {
  # Plots the 2d discrepancy data nicely as different contour plots
  
  library(latex2exp)
  
  # compute mean and stdev of the discrepancy (over the repetitions)
  mean.discr <- apply(discr.data, 1:2, mean)
  stdev.discr <- apply(discr.data, 1:2, sd)
  min.discr <- apply(discr.data, 1:2, min)
  
  # find and plot the minimiser
  min.mean.discr = min(mean.discr)
  print(min.mean.discr)
  min.mean.ind <- arrayInd(which.min(mean.discr), dim(mean.discr))
  print(min.mean.ind)
  
  # find the minimiser of the realised discrepancies
  min.min.discr = min(min.discr)
  print(min.min.discr)
  min.discr.ind <- arrayInd(which.min(min.discr), dim(min.discr))
  print(min.discr.ind)
  
  if (length(discr.data) == 1) {
    warning('Discrepancy data contains only one dimension, not plotting it.')
    
  } else {
    #pdf(file=figure.filename, width = 4, height = 3.5, pointsize = 10, bg = "white")
    pdf(file=figure.filename)
    par(mfrow=c(1,2)) 
    
    # 1/2: plot the mean of the discrepancy vs. neff and mu
    require(grDevices) # for colours
    titl1 <- paste('mean of the "', discrepancy.type, '" discrepancy', sep = '')
    filled.contour(x=neff.grid, y=mu.grid, z=mean.discr, nlevels = 50, color = terrain.colors, 
                   plot.title = title(titl1, xlab=TeX('$n_eff$'), ylab=TeX('$\\mu$')),
                   plot.axes={axis(1); axis(2); points(neff.grid[min.mean.ind[1]],mu.grid[min.mean.ind[2]])})
    
    #lines(x=neff.grid[min.mean.ind[1]], y=mu.grid[min.mean.ind[2]], type = 'p')
    
    # 2/2: plot the stdev of the discrepancy vs. neff and mu
    titl2 <- paste('stdev of the "', discrepancy.type, '" discrepancy', sep = '')
    filled.contour(x=neff.grid, y=mu.grid, z=stdev.discr, nlevels = 50, color = terrain.colors, 
                   plot.title = title(titl2, xlab=TeX('$n_eff$'), ylab=TeX('$\\mu$')))
    
    dev.off()
    if (print.output) {
      print('Discrepancy figure saved.')
    }
  }
  
  # plot raw data values
  if (print.output) {
    print('Mean of the discrepancy:')
    print(mean.discr)
  }
  return(mean.discr)
}


return.abc.post.samples <- function(discr.data, opt2) {
  # Computes the ABC posterior using the table of parameter-discrepancy values
  # In other words, run the rejection sampler ABC algorithm using the precomputed values. 
  
  # discr.data contains now the discrepancy values [neff,mu,rep] -> discrepancy
  dims <- dim(discr.data)
  abc.post <- list()
  abc.post$samples <- matrix(0,dims[1],dims[2]) # each element tells how many simulations are <= epsilon
  
  # determine threshold epsilon, different strategies to choose from...
  if (opt2$epsilon.strategy == 'quantile') {
    epsilon <- quantile(discr.data[is.finite(discr.data)], probs = opt2$epsilon.quantile)
    
  } else if (opt2$epsilon.strategy == 'nr.samples<eps') {
    # choose epsilon so that 'opt2$abc.nr.samples' values are <= epsilon
    sorted.discr <- sort(discr.data[is.finite(discr.data)])
    # check just in case:
    if (length(sorted.discr) < opt2$abc.nr.samples) {
      stop('More abc samples requested than simulations.')
    }
    epsilon <- sorted.discr[opt2$abc.nr.samples]
    
  } else if (opt2$epsilon.strategy == 'fixed') {
    # epsilon is fixed here -> nothing to do
    epsilon <- opt2$epsilon.fixed
    
  } else {
    stop('Incorrect ABC threshold strategy.')
  }
  
  # choose the samples
  abc.post$samples <- apply(discr.data, 1:2, function(x){sum(x<=epsilon)})
  abc.post$epsilon <- epsilon
  
  # NOTE: If the epsilon is set too low, then (in the extreme case) the sample set is empty!
  if (sum(abc.post$samples) <= 0) {
    stop('Sample set is empty.')
  }

  return(abc.post)
}


plot.abc.table <- function(abc.sample.data, neff.grid, mu.grid, titl, figure.filename, 
                           print.output = FALSE, plot.filled.contour = TRUE, lvls = NULL) {
  # Plots the 2d grid-based ABC posterior density 
  
  if (length(abc.sample.data) == 1) {
    warning('ABC data contains only one dimension, not plotting it.')
  } else {
    
    library(latex2exp)
    
    # normalize
    abc.sample.data <- abc.sample.data / sum(abc.sample.data)
    
    # save to file if filename is given
    if (length(figure.filename)>0) {
      #pdf(file=figure.filename, width = 4, height = 3.5, pointsize = 10, bg = "white")
      pdf(file=figure.filename)
    }

    if (plot.filled.contour) {
      
      #library(grDevices) # for contour colours
      #rbPal <- colorRampPalette(c('white','blue'))
      rbPal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      if (is.null(lvls)) {
        nb <- 40
        m <- max(abc.sample.data)
        #print(sum(abc.sample.data)) # should be 1
        lvls <- seq(0,m,len=nb)
      } else {
        nb <- length(lvls)
        # otherwise use the provided levels (needed when different ABC posteriors are compared to ensure
        # that the color lines match!)
      }
      cols <- rbPal(nb)[as.numeric(cut(lvls,breaks = nb))]
      
      # uses a custom contour plot
      filled.contour3(x=neff.grid, y=mu.grid, z=abc.sample.data, #nlevels = 20,  
                      levels = lvls, col = cols, 
                      plot.title = title(titl, xlab=TeX('$n_{eff}$'), ylab=TeX('$\\mu$')), 
                      plot.axes = {axis(side=1,tck = -0.025); axis(side=2,tck = -0.025)})
    } else {
      
      # plot a simple contour plot
      if (is.null(lvls)) {
        nb <- 15
        m <- max(abc.sample.data)
        #print(sum(abc.sample.data)) # should be 1
        lvls <- seq(0,m,len=nb)
      } else {
        nb <- length(lvls)
        # otherwise use the provided levels (needed when different ABC posteriors are compared to ensure
        # that the color lines match!)
      }
      #rbPal <- colorRampPalette(c('white','blue'))
      rbPal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      
      cols <- rbPal(nb)[as.numeric(cut(lvls,breaks = nb))]
      
      contour(x = neff.grid, y = mu.grid, z = abc.sample.data, xaxs = "i", yaxs = "i", tck = -0.025, 
              levels = lvls, drawlabels = FALSE, col = cols)
      title(titl, xlab=TeX('$n_{eff}$'), ylab=TeX('$\\mu$'))
      #abline(v = neff.grid, h = mu.grid)
      #grid(neff.grid, mu.grid)
      grid()
    }
    
    if (length(figure.filename)>0) {
      dev.off()
    }
    if (print.output) {
      print('ABC figure saved.')
    }
  }
  
  # plot raw data values
  if (print.output) {
    print('ABC sample data:')
    print(abc.sample.data)
  }
}

#########################################################################################
#########################################################################################

prior.neff.mu.value <- function(neff, mu, neff.grid, mu.grid, plot.area=FALSE, check.grid=TRUE, root) {
  # Returns 1 if parameter values are inside the prior bounds, and 0 otherwise.
  # This function is used to avoid running simulations with bad parameters (e.g.
  # parameters that are known to cause large discrepancies while requiring huge 
  # computation time.)
  #
  # Check that n_eff and mu are such that they do not produce too many mutations
  # and have thus zero posterior value (and can be directly ignored without 
  # running the costly simulations).
  # This area is based on preliminary simulations and manual inspections. 
  
  #which.curve <- 'poly' # 'poly' or '1/x'
  which.curve  <- '1/x'
  
  points.neff <- c(500,1000,2000,4000,6000,8000,10000) 
  points.mu <- c(0.007,0.004,0.0025,0.0015,0.001,rep(0.0008,2))
  
  if (which.curve == 'poly') {
    points.neff <- points.neff / 10000
    #A <- cbind(points.mu^2, points.mu, 1)
    A <- cbind(points.neff^2, points.neff, 1)
    coeff <- solve(t(A) %*% A, t(A) %*% points.mu) # assuming non-singularity
    #coeff <- solve(A, points.mu)
    coeff = coeff / c(10000^2,10000,1)
    points.neff <- 10000*points.neff
    
  } else {
    A <- cbind(1/points.neff, 1)
    coeff <- solve(t(A) %*% A, t(A) %*% points.mu) # assuming non-singularity
  }
  #print(coeff)
  #if (which.curve  == 'poly' && coeff[1] <= 0) {
  #  stop('Curve fitting for prior density went wrong.')
  #}
  
  if (plot.area) {
    #mu.plot <- seq(min(mu.grid),max(mu.grid),len=100)
    #neff.plot <- coeff %*% rbind(mu.plot^2, mu.plot, 1)
    
    library(latex2exp)
    graphics.off()
    
    neff.plot <- seq(min(neff.grid),max(neff.grid),len=100)
    if (which.curve == 'poly') {
      mu.plot <- t(coeff) %*% rbind(neff.plot^2, neff.plot, 1)
    } else {
      mu.plot <- t(coeff) %*% rbind(1/neff.plot, 1)
    }
    
    figure.filename <- paste(root,'/figures_grids/',sep='')
    figure.filename <- paste(figure.filename,'neff_mu_grid_',length(neff.grid),'x',length(mu.grid),
                             '_',substr(which.curve,1,1),'.pdf',sep='')
    #dev.new()
    pdf(file=figure.filename)
    
    plot(neff.plot, mu.plot, type = 'l', xlab = TeX('$n_{eff}$'), ylab = TeX('$\\mu$'),
         xlim = c(min(neff.grid),max(neff.grid)), ylim = c(min(mu.grid),max(mu.grid)))
    grid()
    lines(neff, mu, type = 'p', pch = 19)
    lines(points.neff, points.mu, type = 'p')
    
    # plot those grid points that are inside bounds
    nr.valid.grid.pts <- 0
    for (i in 1:length(neff.grid)) {
      for (j in 1:length(mu.grid)) {
        if (which.curve == 'poly') {
          #bound.cond <- t(coeff) %*% c(mu^2,mu,1) <= neff
          bound.cond <- t(coeff) %*% c(neff.grid[i]^2, neff.grid[i], 1) <= mu.grid[j]
        } else {
          #bound.cond <- t(coeff) %*% c(1/mu,1) <= neff
          bound.cond <- t(coeff) %*% c(1/neff.grid[i], 1) <= mu.grid[j]
        }
        if (!bound.cond) {
          points(x = neff.grid[i], y = mu.grid[j], pch = 23)
          nr.valid.grid.pts <- nr.valid.grid.pts + 1
        }
      }
    }
    if (TRUE) {
      title(main = paste('Number of grid points = ', nr.valid.grid.pts, sep=''))
    }
    print(paste('Full grid size = ', length(neff.grid)*length(mu.grid), sep=''))
    print(paste('Number of valid grid points = ', nr.valid.grid.pts, sep=''))
    dev.off()
  }
  
  # Check that n_eff and mu are found in the grid.
  if (check.grid && (length(which(mu==mu.grid)) < 1 || length(which(neff==neff.grid)) < 1)) {
    warning('Not a grid point.')
    return(0)
  } else if (!check.grid && (neff < neff.grid[1] || neff > neff.grid[length(neff.grid)] || 
      mu < mu.grid[1] || mu > mu.grid[length(mu.grid)])) {
    return(0)
  }
  # Check that n_eff and mu are inside the prior bounds.
  if (which.curve == 'poly') {
    #bound.cond <- t(coeff) %*% c(mu^2,mu,1) <= neff
    bound.cond <- t(coeff) %*% c(neff^2, neff, 1) <= mu
  } else {
    #bound.cond <- t(coeff) %*% c(1/mu,1) <= neff
    bound.cond <- t(coeff) %*% c(1/neff, 1) <= mu
  }
  if (bound.cond) {
    return(0)
  }
  return(1)
} 




