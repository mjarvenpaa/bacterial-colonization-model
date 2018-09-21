# These functions are used for fitting the mixture model to real or simulated data.

compute.true.data <- function(root) {
  # Fits the mixture model to the MRSA data and plots the figures for the article. 
  
  code.path <- paste(root, '/code', sep='')
  source(paste(code.path, '/mixture_model.R', sep=''))
  source(paste(code.path, '/fit_mixture_model.R', sep=''))
  source(paste(code.path, '/mixture_model_diags.R', sep=''))
  
  simulation.names <- list('mixture_true_data_E_final_pr','mixture_true_data_E_final')
  
  res.dir.name <- paste(root, '/simulation_outputs_mixture_model/', sep='')
  
  save.output.printing <- F # true if terminal output logged to file for collecting the final results
  plot.pred <- T # true if plot the results etc.
  plot.strain.all <- T # true if plot e.g. E and D to same figure for the paper
  
  res <- list()
  i <- 1
  for (simulation.name in simulation.names) {
    if (save.output.printing) {
      con <- file(paste(res.dir.name, simulation.name, '.log', sep = ''))
      sink(con, append = TRUE)
      sink(con, append = TRUE, type = "message")
    }
    
    # fit the mixture model
    fit.mixture.model(root, simulation.name)
    
    # run mcmc diagnostics, compute estimates, run posterior predictive checks etc.
    res[[i]] <- analyse.mixt.model.samples(root, simulation.name, plot.res=plot.pred, plot.pred=plot.pred)
    
    if (save.output.printing) {
      sink()
      sink(type = "message")
    }
    i <- i + 1
  }
  
  # plot all same strain probs for a new observation to one figure
  if (plot.strain.all) {
    filename <- paste(res.dir.name, '/same_strain_prob_all.pdf', sep='')
    plot.strain.pred.figures(res, filename)
  }
}

#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

compute.demos <- function(root) {
  # Fits the mixture model to the artificial data sets and plots the figures for the article. 
  
  code.path <- paste(root, '/code', sep='')
  source(paste(code.path, '/mixture_model.R', sep=''))
  source(paste(code.path, '/fit_mixture_model.R', sep=''))
  source(paste(code.path, '/mixture_model_diags.R', sep=''))
  
  
  ## "N-demo":
  ############
  run.N.demo <- T
  if (run.N.demo) {
    N <- c(50,200,500)
    reps <- 2
    omega1 <- 0.8
    # other parameters are specified elsewhere
    
    for (i in 1:reps) {
      neffs <- list()
      mus <- list()
      titles <- list()
      true.params <- list()
      for (j in 1:length(N)) {
        simulation.name.ij <- paste('artificial_demo_N=',N[j],'_w=',100*omega1,'_rep=',i,sep='')
        
        # fit the mixture model
        fit.mixture.model(root, simulation.name.ij, fake.data.seed=i)
        
        # run mcmc diagnostics, compute estimates, run posterior predictive checks etc.
        res.ij <- analyse.mixt.model.samples(root, simulation.name.ij, plot.res=0, plot.pred=FALSE)
        neffs[[j]] <- res.ij$neffs
        mus[[j]] <- res.ij$mus
        titles[[j]] <- paste('N=', N[j], sep = '')
        true.params[[j]] <- res.ij$clear.data$true.params
      }
      
      # plot the n_eff and mu samples to the same figure (to the last simulation place)
      res.dir.name <- paste(root, '/simulation_outputs_mixture_model/', simulation.name.ij, sep='')
      if (0) {
        filename <- paste(res.dir.name, '/neff_mu_samples_all.pdf', sep='')
        plot.neff.mu.joint.post(neffs, mus, res.ij$neff.grid.new, res.ij$mu.grid.new, plot.grid = TRUE, fig.name = filename, 
                                plot.density.kde = TRUE, plot.first.point = FALSE, true.params = res.ij$clear.data$true.params, 
                                minimal.grid = FALSE, titles = titles)
        
        # plot only in minimal grid
        filename.min <- paste(res.dir.name, '/neff_mu_samples_all_minimal.pdf', sep='')
        plot.neff.mu.joint.post(neffs, mus, res.ij$neff.grid.new, res.ij$mu.grid.new, plot.grid = TRUE, fig.name = filename.min, 
                                plot.density.kde = TRUE, plot.first.point = FALSE, true.params = true.params[[1]], 
                                minimal.grid = TRUE, titles = titles)
      }
    }
  }
  
  
  ## "omega-demo":
  ################
  run.omega.demo <- T
  if (run.omega.demo) { 
    N <- 150
    reps <- 3
    #true.omega1 <- c(0.1,0.25,0.5,0.75,0.9)
    true.omega1 <- c(0.05,0.25,0.5,0.75,0.95)
    # other parameters are specified elsewhere
    
    true.omega1 <- matrix(rep(true.omega1,reps),reps,byrow = TRUE)
    estim.omega1 <- matrix(NA, nrow(true.omega1), ncol(true.omega1))
    true.z <- matrix(NA,N,length(true.omega1))
    estim.z <- matrix(NA,N,length(true.omega1))
    k <- 1
    for (i in 1:ncol(true.omega1)) {
      for (j in 1:reps) {
        simulation.name.ij <- paste('artificial_demo_N=',N,'_w=',100*true.omega1[j,i],'_rep=',j,sep='')
        
        # fit the mixture model
        fit.mixture.model(root, simulation.name.ij, fake.data.seed=j)
        
        # run mcmc diagnostics, compute estimates, run posterior predictive checks etc.
        res.ij <- analyse.mixt.model.samples(root, simulation.name.ij, plot.res=0, plot.pred=FALSE)
        
        estim.omega1[j,i] <- res.ij$omega1$statistics[['Mean']] # posterior mean for omega1
        
        true.z[,k] <- res.ij$clear.data$true.params$z1
        estim.z[,k] <- res.ij$z[,1] # posterior mean for z
        k <- k + 1
      }
    }
    
    # plot the results i.e. omega comparison (to the last simulation place)
    res.dir.name <- paste(root, '/simulation_outputs_mixture_model/', simulation.name.ij, sep='')
    if (0) {
      filename.w <- paste(res.dir.name, '/omega_comparison.pdf', sep='')
      plot_omega_comparison(true.omega1, estim.omega1, filename.w)
    }
  }
  
  
  if (1 && run.N.demo && run.omega.demo) {
    # Plot both of above once more to the same figure
    res.dir.name <- paste(root, '/simulation_outputs_mixture_model/', sep='')
    filename.comb <- paste(res.dir.name, '/N_and_omega.pdf', sep='')
    plot.N.omega.figure(neffs, mus, res.ij$neff.grid.new, res.ij$mu.grid.new, true.params, titles, 
                        true.omega1, estim.omega1, filename = filename.comb, two.row.fig = FALSE)
    filename.comb2 <- paste(res.dir.name, '/N_and_omega2.pdf', sep='')
    plot.N.omega.figure(neffs, mus, res.ij$neff.grid.new, res.ij$mu.grid.new, true.params, titles, 
                        true.omega1, estim.omega1, filename = filename.comb2, two.row.fig = TRUE)
  }
}


#################################################################################################

plot.N.omega.figure <- function(neffs, mus, neff.grid.new, mu.grid.new, true.params, titles, 
                                true.omega1, estim.omega1, filename, two.row.fig = FALSE) {
  # Wrapper for plotting the 'fake' data illustrations i.e. N=50,N=200 plots etc. + omega plot into one figure
  # either as 1x4 or 2x2 plot
  
  # set either to be 1x4 or 2x2 figure if 3+1 plots etc.
  n.plots <- length(neffs) + 1
  if (!two.row.fig) {
    # one row
    pdf(file = filename, width = n.plots*4, height = 4, pointsize = 12)
    par(mfrow = c(1, n.plots)) # 1x...
    par(mai = c(0.6, 0.58, 0.3, 0.2))
  } else {
    # two rows
    pdf(file = filename, width = ceiling(n.plots/2)*4, height = 2*4, pointsize = 12)
    par(mfrow = c(2, ceiling(n.plots/2))) # 2x...
    par(mai = c(0.75, 0.7, 0.3, 0.2))
  }
  
  # Plot N-plot, do not give filename as input so it just plots but does not save the figure
  plot.neff.mu.joint.post(neffs, mus, neff.grid.new, mu.grid.new, plot.grid = TRUE, fig.name = NULL, 
                          plot.density.kde = TRUE, plot.first.point = FALSE, true.params = true.params[[1]], 
                          minimal.grid = TRUE, titles = titles)
  
  # Plot omega-plot to final place
  plot_omega_comparison(true.omega1, estim.omega1, filename=NULL)
  dev.off()
  invisible(0)
}


plot.strain.pred.figures <- function(res, filename) {
  # Wrapper for plotting the same strain posterior probabilities (under different data sets and priors) 
  # for a new measurement (t^*,d^*)
  
  n.plots <- length(res)
  pdf(file = filename, width = n.plots*4, height = 4, pointsize = 12)
  par(mfrow = c(1, n.plots)) 
  if (n.plots == 1) {
    par(mai = c(0.6, 0.6, 0.2, 0.2), mgp = c(1.7,0.55,0))
  } else {
    par(mai = c(0.8, 0.8, 0.3, 0.2), mgp = c(1.7,0.55,0))
  }
  
  # go through different experiments
  s <- 5000
  as <- c('A','B','C','D','E','F')
  for (i in 1:n.plots) {
    plot.strain.pred(res[[i]]$neffs, res[[i]]$mus, res[[i]]$lambdas, res[[i]]$omegas, res[[i]]$zs, res[[i]]$t_0s, 
                     res[[i]]$neff.grid.new, res[[i]]$mu.grid.new, res[[i]]$p.sim, res[[i]]$kt, 
                     d_is = res[[i]]$clear.data$d_is, t_is = res[[i]]$clear.data$t_is, 
                     s = s, filename = NULL, plot.filled = FALSE)
    if (n.plots > 1) {
      title(main = as[i])
    }
  }
  dev.off()
  invisible(0)
}


#################################################################################################

plot_omega_comparison <- function(true.omega1, estim.omega1, filename=NULL) {
  # Plots a comparison of true vs. estimated omega for illustrating that the estimation 
  # algorithm itself works.
  
  library(latex2exp)
  true.omega1 <- as.vector(true.omega1)
  estim.omega1 <- as.vector(estim.omega1)
  
  # draw the figure of those
  if (length(filename) >= 1) {
    le <- 4
    pdf(file=filename, width = le, height = le, pointsize = 12)
    par(mai=c(0.9,0.9,0.2,0.2))
  }
  plot(true.omega1, estim.omega1, type = 'p', pch = 21, col = 'blue',
       xlim = c(0,1), ylim = c(0,1), xaxs="i", yaxs="i", 
       xlab = TeX('$\\omega_{S, true}$'), ylab = TeX('$\\omega_{S, estim}$'), cex.lab = 1.3, main = '')
  lines(c(0,1), c(0,1), type = 'l', col = 'grey')
  grid()
  if (length(filename) >= 1) {
    dev.off()
  }
  invisible(0)
}





