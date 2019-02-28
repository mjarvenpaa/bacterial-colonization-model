# Functions related to the computation of the conditional densities p_S(d_1i|n_eff,mu) and some related computations.

get.p.sim.given.neff.mu.densities <- function(root, opt2, neff.grid, mu.grid, figure.filename = NULL, 
                                              figure.filename2 = NULL, figure.filename3 = NULL, gen.artificial.dens = F) {
  # Loads the raw simulations for computing the conditional densities p_S(d_1i|n_eff,mu) and returns these (normalised)
  # densities.
  
  use.caching <- TRUE # loads a file with precomputed values instead of original raw files
  
  if (gen.artificial.dens) {
    # Generates artificial conditional densities, i.e. density is 0 with probability one
    # This is intended for testing the sampling algorithm only!
    p.sim.given.neff.mu <- rep(list(rep(list(1),opt2$n2)),opt2$n1)
    warning('Using artificial conditional densities of d_i1.', immediate. = TRUE)
    
  } else {
    
    cond.dens.dir.name <- paste(root, '/simulation_outputs_cond_dens/', opt2$simulation.name, sep='')
    precomp.filename <- paste(cond.dens.dir.name, '/precomp_density.RData', sep = '')
    
    # check if the cond.dens have already been computed and cached, otherwise compute it now
    if (use.caching && file.exists(precomp.filename)) {
      
      load(file = precomp.filename) # loads 'di1.data' and 'p.sim.given.neff.mu'
      
    } else {
      
      print('Loading simulations and computing the p_S density...')
      
      ids <- 1:opt2$parallel.runs
      di1.data <- array(NA, c(opt2$n1, opt2$n2, opt2$nr.samples))
      
      for (id.i in ids) {
        cond.dens.filename.id <- paste(cond.dens.dir.name, '/cond_dens_id', id.i, '.RData', sep='')
        load(file = cond.dens.filename.id) # d_i1.samples, neff.grid.cd, mu.grid.cd
        reps <- get.current.ids(id.i, opt2$parallel.runs, opt2$nr.samples)
        if (length(dim(d_i1.samples)) == 3) {
          # in our first simulations, we computed the distances only at one time point (when stability reached)
          di1.data[,,reps] <- d_i1.samples[,,1:length(reps)]
        } else {
          # DISTANCE SAVED AT MULTIPLE TIME POINTS, NOW WE JUST TAKE THE LAST...
          #simul.time.ind <- 1
          simul.time.ind <- dim(d_i1.samples)[4]
          di1.data[,,reps] <- d_i1.samples[,,1:length(reps),simul.time.ind]
        }
      }
      
      ## Go through each 2d grid point and normalise the densities
      # TODO: COULD BE MADE FASTER, PERHAPS
      # NOTE: TIME IS ASSUMED FIXED HERE I.E. ALL THE D_I1 VALUES ARE FROM THE SAME TIME WHEN THE STABILITY IS ACHIEVED
      max.di1 <- max(di1.data)
      p.sim.given.neff.mu <- rep(list(list()),opt2$n1)
      for (i in 1:dim(di1.data)[1]) {
        for (j in 1:dim(di1.data)[2]) {
          if (any(!is.finite(di1.data[i,j,]))) {
            # Simulations were not done because we were outside prior support set
            p.sim.given.neff.mu[[i]][[j]] <- NA
            next
          }
          di1.bins <- -1:max(di1.data[i,j,]) # bins are of the form (a,b] and the first bin thus contains 0's
          h <- hist(di1.data[i,j,], breaks = di1.bins, plot = FALSE)
          di1.counts <- h$counts
          p.sim.given.neff.mu[[i]][[j]] <- di1.counts / sum(di1.counts) # must now be a valid pmf
        }
      }
      
      # cache the computed densities
      if (use.caching) {
        save(di1.data, p.sim.given.neff.mu, file = precomp.filename)
      }
      
      print('Done.')
    }
    
    ## Plot and save the densities p(d_i1|mu,n_eff)
    if (length(figure.filename) > 0) {
      p.max.inds <- 100 # upper bound of plots to draw 
      p.max.inds <- floor(sqrt(p.max.inds))
      
      p.neff.inds <- unique(floor(seq(1,opt2$n1,len=p.max.inds)))
      p.mu.inds <- unique(floor(seq(1,opt2$n2,len=p.max.inds)))
      plot.cond.dens.table(di1.data[p.neff.inds,p.mu.inds,], neff.grid[p.neff.inds], mu.grid[p.mu.inds], 
                           figure.filename = figure.filename, print.info = TRUE)
    }
    
    ## Plot and save another plot with only few densities plotted
    if (length(figure.filename2) > 0) {
      plot.nice.di1.figure(di1.data, neff.grid, mu.grid, figure.filename2)
    }
    
    
    ## Plot also the mean of d_i1 values as a function of n_eff and mu
    if (length(figure.filename3) > 0) {
      plot.di1.means.table(di1.data, neff.grid, mu.grid, figure.filename3)
    }
  }
  return(p.sim.given.neff.mu)
}


plot.cond.dens.table <- function(d_i1.samples, neff.grid, mu.grid, figure.filename, print.info=TRUE) {
  # Plots the simulated densities (separately for each grid point) and saves them to file
  # That is, plots the estimated densities p_S(d_1i|n_eff,mu) for different n_eff,mu values
  
  n.neff.grid <- length(neff.grid)
  n.mu.grid <- length(mu.grid)
  
  di1.max <- max(d_i1.samples[is.finite(d_i1.samples)])
  col.plots <- 2
  
  # compute number of plots to draw
  #n.tot <- n.mu.grid * n.neff.grid
  n.tot <- 0
  for (i in 1:n.neff.grid) {
    for (j in 1:n.mu.grid) {
      if (all(is.finite(d_i1.samples[i,j,]))) {
        n.tot <- n.tot + 1
      }
    }
  }
  
  if (n.tot > 200) {
    stop('Trying to draw too many plots.')
  } else if (n.tot <= 0) {
    return()
  }
  
  hei <- ceiling(n.tot/col.plots)
  pdf(file=figure.filename, width = 3*col.plots, height = 2.25*hei)
  par(mfrow=c(hei,col.plots))
  
  for (i in 1:n.neff.grid) {
    for (j in 1:n.mu.grid) {
      if (all(is.finite(d_i1.samples[i,j,]))) {
        # we are inside prior support so draw the plot!
        titl.ij <- paste('p(d_i1|n_eff=', neff.grid[i], ', mu=', signif(mu.grid[j],digits=3), ')', sep = '')
        hist(d_i1.samples[i,j,], main = titl.ij, xlab = 'd_i1', ylab = 'Simulated freq.', breaks = 0:(di1.max+1), right=F,
             col = 'lightblue')
      }
    }
  }
  dev.off()
  
  # if (print.info) {
  #   #print('d_1i counts over all the distances')
  #   #print(table(d_i1.samples))
  #   #print('')
  #   n_eff.ind <- 1
  #   mu.wise.counts <- apply(d_i1.samples[n_eff.ind,,],1,table)
  #   names(mu.wise.counts) <- signif(mu.grid, digits = 5)
  #   print('d_1i counts')
  #   print(mu.wise.counts)
  # }
}


plot.simul.cond.dens <- function(d_i1.samples, neff.param, mu.param, figure.filename) {
  # Plots essentially a similar figure as the above function but in a bit nicer way. 
  
  library(latex2exp)
  
  di1.max <- max(d_i1.samples[is.finite(d_i1.samples)])
  y <- nrow(neff.param)
  x <- ncol(neff.param)
  
  pdf(file=figure.filename, width = 1.2*2.2*x, height = 1.2*1.8*y)
  par(mfrow=c(y,x),tcl=-0.4)
  
  for (i in 1:y) {
    for (j in 1:x) {
      titl.ij <- sprintf('$p(d_{i1}|n_{eff} = %i,\\, \\mu = %f)$', neff.param[i,j], signif(mu.param[i,j],digits=3))
      par(mai=c(0.275,0.01,0.2,0.05))
      if (j == 1) {
        par(mai=c(0.275,0.01,0.2,0.05))
      }
      hist(d_i1.samples[i,j,], main = TeX(titl.ij), xlab = '', ylab = '', breaks = 0:(di1.max+1), right=F, freq=F, yaxt='n',
           col = 'lightblue')
    }
  }
  dev.off()
}


plot.nice.di1.figure <- function(d_i1.samples, neff.grid, mu.grid, figure.filename) {
  # Wrapper for plotting a nice d_i1 density plot figure for the paper.
  
  # select which n_eff,mu-values to plot
  n <- length(neff.grid)
  m <- length(mu.grid)
  x <- 3
  y <- 2
  if (n <= 2 || m <= 2) {
    # testing case, do nothing
    invisible(0)
    
  } else if (n == 40 && m == 40) {
    # our original grid used for computing preliminary results
    #neff.ind <- matrix(c(1,5,n,1,5,n),y,x,byrow = T)
    #mu.ind <- matrix(c(5,5,5,m,m,m),y,x,byrow = T)
    neff.ind <- matrix(c(2,2,10,10,20,n),y,x,byrow = T)
    mu.ind <- matrix(c(2,20,2,10,2,2),y,x,byrow = T)
    
  } else if (n == 100 && m == 100) {
    # full grid for final results
    neff.ind <- matrix(c(5,5,25,25,50,n),y,x,byrow = T)
    mu.ind <- matrix(c(5,50,5,25,10,10),y,x,byrow = T)
    
  } else {
    browser()
  }
  
  neff.param <- matrix(neff.grid[neff.ind],y,x,byrow = F)
  mu.param <- matrix(mu.grid[mu.ind],y,x,byrow = F)
  samples <- array(NA,c(y,x,dim(d_i1.samples)[3]))
  for (i in 1:y) {
    for (j in 1:x) {
      samples[i,j,] <- d_i1.samples[neff.ind[i,j],mu.ind[i,j],]
    }
  }
  plot.simul.cond.dens(samples, neff.param, mu.param, figure.filename)
}


plot.di1.means.table <- function(d_i1.samples, neff.grid, mu.grid, figure.filename) {
  # Plots the mean of d_i1 values as a function of n_eff and mu
  
  library(latex2exp)
  
  # Compute the means first
  mean.d <- rowMeans(d_i1.samples, dims = 2)
  lvls <- seq(0,min(20,ceiling(max(mean.d,na.rm=TRUE))),by=0.5)
  
  pdf(file = figure.filename, width = 6, height = 6)
  contour(x = neff.grid, y = mu.grid, z = mean.d, levels = lvls, xaxs="i", yaxs="i", 
          xlim = c(min(neff.grid),max(neff.grid)), ylim = c(min(mu.grid),max(mu.grid)))
  title(xlab = TeX('$n_{eff}$'), ylab = TeX('$\\mu$'), main = TeX('estimated mean of $d_{i1}$'))
  grid()
  dev.off()
} 




