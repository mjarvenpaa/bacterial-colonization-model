# Analyse the mixture model fitting results (assess whether MCMC chains have converged, compute the final estimates etc.)

analyse.mixt.model.samples <- function(root, simulation.name=NULL, plot.res=TRUE, plot.pred=TRUE) {
  
  if (is.null(simulation.name)) {
    simulation.name <- 'mixture_test_artificial' # USE ARTIFICIAL DATA
  }
  
  ## Settings etc.
  code.path <- paste(root, '/code', sep='')
  source(paste(code.path, '/get_settings.R', sep=''))
  source(paste(code.path, '/get_data.R', sep=''))
  source(paste(code.path, '/abc_methods_simul_model.R', sep=''))
  source(paste(code.path, '/cond_dens_help_functions.R', sep=''))
  source(paste(code.path, '/mixture_model.R', sep=''))
  source(paste(code.path, '/figure_functions.R', sep=''))
  source(paste(code.path, '/auxiliary_functions.R', sep=''))
  
  res.dir.name <- paste(root, '/simulation_outputs_mixture_model/', simulation.name, sep='')
  
  # get settings used in Gibbs sampling
  opt3 <- get.mixture.fitting.settings(simulation.name, print.settings=FALSE)
  kt <- opt3$kt
  nr.chains <- opt3$nr.chains
  nr.gibbs.samples <- opt3$nr.gibbs.samples
  
  opt2 <- get.cond.dens.settings(simulation.name=opt3$cond.dens.simulation.name, print.settings=FALSE)
  neff.grid.new <- opt2$neff.grid
  mu.grid.new <- opt2$mu.grid
  
  
  ## CONVERGENCE CHECKINGS
  
  graphics.off()
  library('coda')
  
  ## Collect samples in the chains
  chains4coda <- list()
  for (cur.chain in 1:nr.chains) {
    
    # Load the current samples
    load(file = paste(res.dir.name, '/mixture_fitting_test', cur.chain, '.RData', sep=''))
    
    # TODO: this is a bit ugly...
    if (cur.chain == 1) {
      # Initialize
      N <- dim(zs)[1]
      final.inds <- ceiling((nr.gibbs.samples+1)/2+1):(nr.gibbs.samples+1) # the first half ignored as burn-in
      chain.size <- length(final.inds)
      final.nr.samples <- nr.chains * chain.size
      
      zs.res <- array(NA,c(N,2,final.nr.samples)) 
      omegas.res <- matrix(NA,2,final.nr.samples) 
      lambdas.res <- rep(NA,final.nr.samples) 
      neffs.res <- rep(NA,final.nr.samples) 
      mus.res <- rep(NA,final.nr.samples) 
      t_0s.res <- matrix(NA,N,final.nr.samples) 
      log.unnorm.post.res <- rep(NA,final.nr.samples)
    }
    z.plot.inds <- get.z.plot.inds(d_is) # which z and t_0 components are plotted
    
    # Combine chains
    inds <- ((cur.chain-1)*chain.size+1):(cur.chain*chain.size)
    zs.res[,,inds] <- zs[,,final.inds]
    omegas.res[,inds] <- omegas[,final.inds] 
    lambdas.res[inds] <- lambdas[final.inds]
    neffs.res[inds] <- neffs[final.inds]
    mus.res[inds] <- mus[final.inds]
    t_0s.res[,inds] <- t_0s[,final.inds]
    log.unnorm.post.res[inds] <- log.unnorm.post[final.inds]
    
    # Thinning (if necessary)
    #...
    
    # Reformat chain for coda
    # NOTE: neff, mu, z not included if just one value e.g. during the preliminary grid computations
    if (length(unique(neffs[final.inds])) < 2 || length(unique(mus[final.inds])) < 2) {
      theta.res <- cbind(lambdas[final.inds], as.double(omegas[1,final.inds]), t(t_0s[,final.inds]))
      warning('Excluded mu and n_eff from convergence checkings.', immediate. = TRUE)
    } else {
      theta.res <- cbind(neffs[final.inds],mus[final.inds],lambdas[final.inds], 
                         as.double(omegas[1,final.inds]), t(t_0s[,final.inds]))
    }
    chains4coda[[cur.chain]] <- as.mcmc(theta.res)
  }

  
  ## Check the convergence of the chains using coda package
  if (nr.chains > 1) {
    chain4coda <- as.mcmc.list(chains4coda)
    cat('Using coda...\n')
    #gelman.plot(chain4coda)
    print(gelman.diag(chain4coda))
    cat('\n')
  }
  
  
  ## Correlation plots, could also just use 'pairs'
  library('BayesianTools') # for correlation plot
  neff <- neffs.res
  mu <- mus.res
  omega_1 <- as.double(omegas.res[1,])
  lambda <- lambdas.res
  t_0 <- t(t_0s.res[z.plot.inds,,drop=FALSE])
  
  sample.data4cor <- cbind(neff, mu, lambda, omega_1, t_0)
  colnames(sample.data4cor) <- c('n_eff','mu','lambda','omega_1',paste('t_0',z.plot.inds,sep=''))
  
  corr.filename <- paste(res.dir.name, '/corr_plot.pdf', sep = '')
  if (plot.res) {
    pdf(file=corr.filename, width=10, height=10, pointsize=12) 
    correlationPlot(sample.data4cor, labels=c('1','2'))
    dev.off()
    
    ## Plot combined chains etc.
    sample.filename.res <- paste(res.dir.name, '/fig_samples_final.pdf', sep='')
    simple.sample.plot(zs.res, omegas.res, lambdas.res, neffs.res, mus.res, t_0s.res, log.unnorm.post.res, 
                       clear.data$t_is, clear.data$d_is, neff.grid.new, mu.grid.new, z.plot.inds, sample.filename.res) 
    
    sample.filename.res2 <- paste(res.dir.name, '/fig_neff_vs_mu_final.pdf', sep='')
    plot.neff.mu.joint.post(neffs.res, mus.res, neff.grid.new, mu.grid.new, fig.name=sample.filename.res2, 
                            plot.first.point=FALSE, true.params=clear.data$true.params)
    sample.filename.res3 <- paste(res.dir.name, '/fig_neff_vs_mu_final_minimal.pdf', sep='')
    plot.neff.mu.joint.post(neffs.res, mus.res, neff.grid.new, mu.grid.new, fig.name=sample.filename.res3, 
                            plot.first.point=FALSE, true.params=clear.data$true.params, minimal.grid = TRUE)
  }
  
  
  ## Print final, estimated values
  qs <- c(2.5,50,97.5)/100
  
  cat('\n')
  cat('N_EFF:')
  summ.neff <- summary(as.mcmc(neffs.res), quantiles=qs)
  print(summ.neff)
  
  cat('MU:')
  summ.mu <- summary(as.mcmc(mus.res), quantiles=qs)
  print(summ.mu)
  
  cat('LAMBDA:')
  summ.lambda <- summary(as.mcmc(lambdas.res), quantiles=qs)
  print(summ.lambda)
  
  cat('OMEGA_1:')
  summ.omega1 <- summary(as.mcmc(omegas.res[1,]), quantiles=qs)
  print(summ.omega1)
  
  ## Z
  z.print.inds <- seq(1, dim(zs.res)[1], by=20)
  z.means <- apply(apply(zs.res, 2:3, as.double), 1:2, mean) # TODO: COULD BE MADE FASTER
  zs.print <- cbind(z.means, t_is, d_is)
  colnames(zs.print) <- c('z1','z2','t_i','d_i')
  rownames(zs.print) <- 1:nrow(zs.print)
  cat('MEAN OF Z:\n')
  print(zs.print[z.print.inds,])
  cat('\n')
  
  ## T_0
  cat('T_0:')
  summ.t0 <- summary(as.mcmc(t(t_0s.res[z.print.inds,])), quantiles=qs)
  print(summ.t0)
  
  
  ## WHAT KIND OF DATA IS PREDICTED
  p.sim <- get.p.sim.given.neff.mu.densities(root, opt2, neff.grid.new, mu.grid.new)
  if (plot.pred) {
    obs.data.filename <- paste(res.dir.name, '/post_pred_data_observed.pdf', sep = '')
    plot.observed.di.density(t_is, d_is, obs.data.filename)
    
    s <- 10000
    t_i.test <- 2000
    filename <- cbind(paste(res.dir.name, '/post_pred_data_same_strain.pdf', sep = ''),
                      paste(res.dir.name, '/post_pred_data_diff_strain.pdf', sep = ''))
    for (i in 1:2) {
      plot.di.density(neffs.res, mus.res, lambdas.res, neff.grid.new, mu.grid.new, p.sim, 
                      kt, t_i = t_i.test, same.strain = (i==1), filename = filename[i], s = s)
    }
  }
  
  ## POSTERIOR PREDICTIVE CHECK
  if (plot.pred) {
    s <- 5000
    filename <- paste(res.dir.name, '/post_pred_vs_data_all_single.pdf', sep = '')
    filename2 <- paste(res.dir.name, '/post_pred_vs_data.pdf', sep = '')
    post.pred.check(neffs.res, mus.res, lambdas.res, zs.res, t_0s.res, neff.grid.new, mu.grid.new, 
                    p.sim, kt, t_is, d_is, filename, filename2, s = s)
  }
  
  # THE HEATMAP TYPE OF ILLUSTRATION OF THE Z-PROB. SURFACE FOR A NEW MEASUREMENT (t^*,d^*)
  if (plot.pred) {
    s <- 5000
    filename <- paste(res.dir.name, '/z_vs_t_and_d.pdf', sep = '')
    plot.strain.pred(neffs.res, mus.res, lambdas.res, omegas.res, zs.res, t_0s.res, neff.grid.new, mu.grid.new, 
                     p.sim, kt, d_is, t_is, s = s, filename = filename)
  }
  
  
  # Return the computed mcmc summaries and combined chains etc.
  invisible(list('clear.data'=clear.data, 'neff'=summ.neff, 'mu'=summ.mu, 'lambda'=summ.lambda, 'omega1'=summ.omega1,
            'z'=z.means, 't0'=summ.t0, 'kt'=kt, 'p.sim'=p.sim, 
            'neffs'=neffs.res, 'mus'=mus.res, 'lambdas'=lambdas.res, 'omegas'=omegas.res, 'zs'=zs.res, 't_0s'= t_0s.res,
            'neff.grid.new'=neff.grid.new, 'mu.grid.new'=mu.grid.new))
}

##########################################################################################################
##########################################################################################################

plot.observed.di.density <- function(t_is, d_is, filename) {
  # Plot the observed d_i-values into one figure.
  
  if (length(filename) > 0) {
    
    library(latex2exp)
    
    titl = 'obs. data, all distances'
    #breaks <- seq(0,(max(d_is)+1),len=min(max(d_is)/5,100))
    breaks <- 0:(max(d_is)+1)
    
    pdf(file = filename)
    hist(x = d_is, xlab = TeX('$d_i$'), main = titl, breaks = breaks, right=F, freq = FALSE)
    dev.off()
  }
}


plot.di.density <- function(neffs, mus, lambdas, neff.grid, mu.grid, p.sim, kt, t_i = 1000, 
                            same.strain = TRUE, filename, s = 1) {
  # Plot the posterior predictive density p(di|t_i) for a given value t_i, based on the fitted model
  # and the given information of whether the strain is the same or not.
  
  library(latex2exp)
  
  # sample posterior parameter(s)
  r.inds <- sample(x = 1:length(neffs), size = s, replace = TRUE)
  if (same.strain) {
    # same strain case, need the simulator model here
    neff <- neffs[r.inds]
    mu <- mus[r.inds]
    d_i2 <- rpois(n = s, lambda = mu * t_i)
    # since we may not have computed the density p(d_i1\n_eff, mu) for the given n_eff, mu value 
    # we need to interpolate here (it should add only small error)
    d_i1 <- rep(NA,s)
    for (i in 1:s) {
      # find the closest grid point at which the density is computed and sample from it
      # the current point used is one of the corner points while one could also interpolate
      ij <- get.2x2.grid(neff.grid, mu.grid, c(neff[i],mu[i]))
      probs <- p.sim[[ij[1]]][[ij[2]]]
      d_i1[i] <- sample(x = 0:(length(probs)-1), size = 1, replace = TRUE, prob = probs)
    }
    d_i <- d_i1 + d_i2
    titl <- 'posterior predictive (same strain case)'
    breaks <- 0:(max(d_i)+1)
    
  } else {
    # different strain case
    lambda <- lambdas[r.inds]
    t_i0 <- rgamma(n = s, shape = kt, rate = lambda)
    mut <- mus[r.inds] * (2*t_i0 + t_i)
    d_i <- rpois(n = s, lambda = mut)
    titl <- 'posterior predictive (different strain case)'
    #breaks <- seq(0,max(d_i)+1,len=min(max(d_i)/5,100))
    breaks <- 0:(max(d_i)+1)
  }
  
  # plot the histogram
  if (length(filename) > 0 && s > 1) {
    
    library(latex2exp)
    
    pdf(file = filename)
    hist(x = d_i, xlab = TeX('$d_i$'), main = titl, breaks = breaks, right=F, freq = FALSE)
    dev.off()
  }
  invisible(d_i)
}


post.pred.check <- function(neffs, mus, lambdas, zs, t_0s, neff.grid, mu.grid, p.sim, kt, 
                            t_is, d_is, filename = NULL, filename2 = NULL, s = 100) {
  # Compute posterior predictive distribution for each data point (i.e. consequtive hospital visits)
  # and plot it wrt. the corresponding data point.
  # Plots also the replicated datasets similar to the observed ones. 
  
  library(latex2exp)
  
  # sample posterior parameter(s)
  #r.inds <- sample(x = 1:length(neffs), size = s, replace = TRUE)
  r.inds <- floor(seq(1,length(neffs),len=s))
  z1 <- zs[,1,r.inds]
  t0 <- t_0s[,r.inds]
  neff <- neffs[r.inds]
  mu <- mus[r.inds]
  
  N <- length(t_is)
  D <- matrix(NA,N,s) # generated distances as N x s matrix!
  for (i in 1:s) {
    model1.inds <- which(z1[,i])
    model2.inds <- which(!z1[,i])
    
    # model 1 case
    # find the closest grid point at which the density is computed and sample from it
    # the current point used is one of the corner points while one could also interpolate
    ij <- get.2x2.grid(neff.grid, mu.grid, c(neff[i],mu[i]))
    probs <- p.sim[[ij[1]]][[ij[2]]]
    di <- sample(x = 0:(length(probs)-1), size = length(model1.inds), replace = TRUE, prob = probs)
    D[model1.inds,i] <- di + rpois(n = length(model1.inds), lambda = mu[i]*t_is[model1.inds])
    
    # model 2 case
    mut <- mu[i] * (2*t0[model2.inds,i] + t_is[model2.inds])
    D[model2.inds,i] <- rpois(n = length(model2.inds), lambda = mut)
  }
  
  ## 1) Plot the predictive density of each point vs. this point
  if (length(filename) > 0) {
    points.to.draw <- seq(1,N,by=3)
    #points.to.draw <- unique(floor(seq(1,N,len=20)))
    
    col.plots <- 3
    hei <- ceiling(length(points.to.draw)/col.plots)
    pdf(file=filename, width = 3*col.plots, height = 2.25*hei)
    par(mfrow=c(hei,col.plots))
    
    par(mai=c(0.275,0.25,0.2,0.001))
    
    #breaks <- 0:(max(d_is,D)+1)
    for (i in points.to.draw) {
      breaks <- 0:(max(D[i,],d_is[i])+1)
      titl <- paste('d_i=', d_is[i], ', t_i=', round(t_is[i]), sep = '')
      hist(x = D[i,], xlab = TeX('$d_i$'), ylab = 'Post. pred. density', main = titl, breaks = breaks, freq = FALSE,
           col = 'lightblue')
      lines(x = d_is[i]*c(1,1), y = c(0,1), col = 'red') # plot also d_i measurement with red
    }
    dev.off()
  }
  
  ## 2) Plot another figure where full observed data is graphically compared to some full simulated replicated data
  if (length(filename2) > 0) {
    reps.to.draw <- 5
    
    ylim.val <- NULL
    ylim.val <- c(0,0.4)
    
    # 2 x floor(N/2)+1 plot: observed data to (1,1), simulated data to other places
    col.plots <- floor(reps.to.draw/2)+1
    row.plots <- 4
    pdf(file=filename2, width = 3*col.plots, height = 2.25*row.plots)
    par(mfrow=c(row.plots,col.plots))
    
    # plot observed data
    par(mai=c(0.275,0.25,0.2,0.001))
    titl = 'observed data'
    breaks <- 0:(max(D,d_is)+1)
    hist(x = d_is, xlab = TeX('$d_i$'), main = titl, breaks = breaks, freq = FALSE,
         col='lightblue', ylim = ylim.val)
    
    # plot replicated data from posterior predictive
    for (i in 1:reps.to.draw) {
      titl.i <- paste('simulated data, rep=',i,sep='')
      hist(x = D[,i], xlab = TeX('$d_i$'), main = titl.i, breaks = breaks, freq = FALSE,
           col='lightblue', ylim = ylim.val)
    }
    
    # 2b) Plot above also as constrained to be <= 50 for better visualisation of small distances
    cut.dist <- 50
    
    #titl = paste('observed data (d_i<',cut.dist,')',sep='')
    titl = 'observed data'
    hist(x = d_is, xlab = TeX('$d_i$'), main = titl, breaks = breaks, freq = FALSE, xlim = c(0,cut.dist),
         col='lightblue', ylim = ylim.val)
    for (i in 1:reps.to.draw) {
      #titl.i <- paste('simulated data (d_i<',cut.dist,'), rep=',i,sep='')
      titl.i <- paste('simulated data, rep=',i,sep='')
      hist(x = D[,i], xlab = TeX('$d_i$'), main = titl.i, breaks = breaks, freq = FALSE, xlim = c(0,cut.dist),
           col='lightblue', ylim = ylim.val)
    }
    dev.off()
  }
  invisible(D)
}


plot.strain.pred <- function(neffs, mus, lambdas, omegas, zs, t_0s, neff.grid, mu.grid, p.sim, kt, d_is = NULL, t_is = NULL, 
                             d.grid = seq(0,50,by=2), t.grid = seq(100,6000,by=200), s = 1000, filename = NULL, plot.filled = FALSE) {
  # Computes and plots the probability of same strain assignment of a new (future) measurement with
  # distance d^* and time interval t^*. These are computed and plotted for a grid of input values.
  # If filename is given then the figure is saved to file. 
  
  # Compute in the grid
  nd <- length(d.grid)
  nt <- length(t.grid)
  n <- length(neffs)
  #r.inds <- sample(x = 1:n, size = s, replace = T) # indexes of the posterior samples
  r.inds <- floor(seq(1,n,len=s))
  t0.stars <- rgamma(n = s, shape = kt, rate = lambdas[r.inds])
  
  z1.star.grid <- array(NA,c(nd,nt,s))
  z2.star.grid <- array(NA,c(nd,nt,s))
  
  dt.all <- meshgrid(d.grid, t.grid)
  d.all.vec <- as.vector(dt.all$x)
  t.all.vec <- as.vector(dt.all$y)
  for (i in 1:s) {
    # model 1 case i.e. z^*_1 == 1
    # find the closest grid point at which the density is computed and evaluate it there
    # the current point used is one of the corner points while one could also interpolate
    ij <- get.2x2.grid(neff.grid, mu.grid, c(neffs[r.inds[i]],mus[r.inds[i]]))
    log.z1.star.grid <- log.p_S.density(d.all.vec, t.all.vec, mus[r.inds[i]], mu.grid, p.sim[[ij[1]]][[ij[2]]])
    z1.star.grid[,,i] <- exp(log(omegas[1,r.inds[i]]) + matrix(log.z1.star.grid, nd, nt))
    
    # model 2 case i.e. z^*_2 == 1
    mut <- mus[r.inds[i]] * (2*t0.stars[i] + t.all.vec)
    z2.star.grid[,,i] <- omegas[2,r.inds[i]] * matrix(dpois(x = d.all.vec, lambda = mut),nd,nt)
  }
  
  # normalise
  # we normalise separately for each case since p(y) is constant only for fixed y!
  z1.star.grid <- 1/s * rowSums((z1.star.grid), dims = 2)
  z2.star.grid <- 1/s * rowSums((z2.star.grid), dims = 2)
  z1.star.grid <- z1.star.grid / (z1.star.grid + z2.star.grid)
  z2.star.grid <- 1 - z1.star.grid
  
  ######### Plot ###########
  
  if (length(filename) > 0) {
    le <- 5
    pdf(file=filename, width = le, height = le, pointsize = 12)
    #par(mai=c(0.275,0.25,0.2,0.001))
  }
  
  # SET UP DATAPOINTS AND THEIR ESTIMATED CLASSES
  jitter <- 0
  if (length(unique(t_is)) < 10) {
    jitter <- 50
  }
  
  y <- t_is + runif(n = length(t_is), min = -jitter, max = jitter)
  rbPal <- colorRampPalette(c('red','blue'))
  #rbPal <- colorRampPalette(c('white','blue'))
  #rbPal <- colorRampPalette(c('darkred','red','blue','lightblue'))
  # compute posterior means of same strain cases
  mean_z1 <- rowMeans(zs[,1,])
  nb <- 20
  cols <- rbPal(nb)[as.numeric(cut(mean_z1,breaks = nb))]

  
  if (!plot.filled) {
    # PLOT THE CONTOUR
    
    # plot datapoints
    plot(NA, xlab = 'd*', ylab = 't* (generations)', 
         xlim = c(d.grid[1]-0.5,d.grid[length(d.grid)]), ylim = c(t.grid[1],t.grid[length(t.grid)]),
         xaxs="i", yaxs="i") 
    lines(x = d_is, y = y, type = 'p', pch = 20, col = cols)
    
    # print also the mean_z
    if (F) {
      print(sort(mean_z1))
    }
    
    # plot the contour plot
    lvls <- c(0.01,0.05,seq(0.1,0.9,by=0.1),0.95,0.99)
    lvls <- c(0.001,lvls,0.999)
    n.lvls <- length(lvls)
    nb <- 20
    cols <- rbPal(nb)[as.numeric(cut(lvls,breaks = nb))]
    contour(x = d.grid, y = t.grid, z = z1.star.grid, levels = lvls, zlim = c(0,1), add = TRUE,
            col = cols)
    grid()
    
  } else {
    # plot filled contour instead
    
    lvls <- seq(0,1,by=0.05)
    n.lvls <- length(lvls)
    filled.contour3(x = d.grid, y = t.grid, z = z1.star.grid, xaxs="i", yaxs="i", levels = lvls,
                      plot.axes={axis(1); axis(2); lines(x = d_is, y = y, type = 'p', pch = 21, 
                                                         bg = cols, col = 'black')},
                      color.palette = rbPal)
    title(xlab = 'd', ylab = 't (generations)')
  }
  
  if (length(filename) > 0) {
    dev.off()
  }
  invisible(z1.star.grid)
}





