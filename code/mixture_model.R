# This file contains some misc functions needed for fitting the mixture model and related plotting

rand.from.prior <- function(N, alpha_lambda, beta_lambda, gamma_omega, neff.grid, mu.grid, kt, log.neff.mu.prior=NULL) {
  # Generates one random initialisation point from the prior of the parameter theta of the mixture model.
  # This is used for initializing the MCMC chains.
  
  r <- as.logical(sample(x = 0:1, size = N, replace = TRUE))
  z <- cbind(r,!r)
  #omega <- c(sum(r),sum(!r)) / length(r)
  # generating from the 2d Dirichlet 
  omega <- rep(0,2)
  omega[1] <- rbeta(n = 1, shape1 = gamma_omega[1], shape2 = gamma_omega[2])
  omega[2] <- 1 - omega[1] 
  lambda <- rgamma(n = 1, shape = alpha_lambda, rate = beta_lambda)
  
  # NOTE: If prior is not given as input, generate uniformly from the grid
  if (all(log.neff.mu.prior==0)) {
    max.nr.tries <- 10000
    for (i in 1:max.nr.tries) {
      neff <- sample(x = neff.grid, size = 1)
      mu <- sample(x = mu.grid, size = 1)
      if (prior.neff.mu.value(neff, mu, neff.grid, mu.grid) == 1) {
        # sampled value is in prior support so take it
        break
      }
      if (i==max.nr.tries) {
        stop('Did not find a prior parameter inside the bounds.') # should not happen with very high prob.
      }
    }
  } else {
    #warning('Generated n_eff, mu uniformly from the grid for the initialisation.', immediate. = TRUE)
    ind <- sample.cat.density.matrix(log.neff.mu.prior)
    neff <- neff.grid[ind[1]]
    mu <- mu.grid[ind[2]]
  }
  
  t_0 <- rgamma(n = N, shape = kt, rate = lambda)
  
  r.prior <- list(z=z, omega=omega, lambda=lambda, neff=neff, mu=mu, t_0=t_0)
  return(r.prior)
}


p_S.density <- function(d_i, t_i, mu, p.sim.given.neff.mu) {
  # Computes the value of p_S(d_i|mu,t_i) for given d_i, mu and t_i.
  # d_i and t_i can be vectors, mu must be constant and p.sim.given.mu is the p_sim(d_i1|mu,data_1:2)
  # density estimated with ABC.
  # NOTE: Currently this is not needed for anything since all computations are done in log-domain.
  
  N <- length(d_i)
  muti <- mu * t_i
  d_max <- length(p.sim.given.mu) - 1 # first index is zero distance
  p <- rep(NA,N)
  
  for (j in 1:N) {
    # uses the convolution formula for the sum of two discrete rv's
    ks.j <- max(d_i[j] - d_max,0):d_i[j]
    #poi.terms <- muti[j]^(ks.j) * exp(-muti[j]) / factorial(ks.j)
    log.poi.terms <- (ks.j) * log(muti[j]) - muti[j] - lfactorial(ks.j)
    p_sim.terms <- p.sim.given.neff.mu[d_i[j] - ks.j + 1]
    # This can cause underflow here:
    #p[j] <- sum(exp(log.poi.terms) * p_sim.terms)
    p[j] <- sum(exp(log.poi.terms + log(p_sim.terms)))
  }
  return(p)
}


log.p_S.density <- function(d_i, t_i, mu, mu.grid, p.sim.given.neff.mu) {
  # Computes the value of log(p_S(d_i|neff,mu,t_i)) for given d_i, n_eff, mu and t_i.
  # d_i and t_i can be vectors, neff and mu must be constants, neff.grid and mu.grid must be grid vectors 
  # containing the corresponding neff and mu values and p.sim.given.mu is the matrix containing the 
  # p_sim(d_i1|n_eff,mu,data_1:2) values. 
  # This log version should produce more accurate and robust computations!
  
  N <- length(d_i)
  muti <- mu * t_i
  d_max <- length(p.sim.given.neff.mu) - 1 # first index is zero distance
  log.p <- rep(NA,N)
  
  for (j in 1:N) {
    # uses the convolution formula for the sum of two discrete rv's
    ks.j <- max(d_i[j] - d_max,0):d_i[j]
    log.poi.terms <- (ks.j) * log(muti[j]) - muti[j] - lfactorial(ks.j)
    p_sim.terms <- p.sim.given.neff.mu[d_i[j] - ks.j + 1]
    log.p[j] <- log.sum.exp(log.poi.terms + log(p_sim.terms))
  }
  return(log.p)
}


log.sum.exp <- function(log.terms) {
  # We use the identity log(sum(exp(x))) = max(x) + log(sum(exp(x - max(x)))) to compute this value
  # so that possible under/overflows are handled more robustly than using the naive implementation.
  # Works for both vector and matrix input. 
  
  #max.log.term <- max(log.terms)
  max.ind <- which.max(log.terms)
  max.ind <- max.ind[1]
  max.log.term <- log.terms[max.ind]
  log.sum <- max.log.term + log1p(sum(exp(log.terms[-max.ind] - max.log.term)))
  return(log.sum)
}


sample.cat.density.matrix <- function(log.probs.matrix) {
  # Draws a sample from a categorical density when logarithm of each unnormalised bin probability
  # is given as input. 
  # Gumbel max trick is used for computations to avoid underflow.
  # This version assumes that log.probs is a matrix and returns the row and column corresponding 
  # to the randomly chosen category.
  
  k <- dim(log.probs.matrix)
  u_is <- matrix(runif(k[1]*k[2]), k[1], k[2])
  g_is <- -log(-log(u_is))
  r <- g_is + log.probs.matrix
  r <- which(r == max(r), arr.ind = TRUE)
  s <- r[1,] # duplicates happen with probability zero...
  return(s) # each element contains the two indexes
}


sample.cat.density <- function(log.probs) {
  # Draws a sample from a categorical density when logarithm of each unnormalised bin probability
  # is given as input. 
  # Gumbel max trick is used for computations to avoid underflow.
  
  k <- length(log.probs)
  u_is <- runif(k)
  g_is <- -log(-log(u_is))
  s <- which.max(g_is + log.probs)
  return(s) # returns a scalar index
}


sample.cat.density.vectorized <- function(log.probs) {
  # Draws samples from a categorical density when logarithm of each unnormalised bin probability
  # is given as input. 
  # Gumbel max trick is used for computations to avoid underflow.
  # Here each row in log.probs is computed separately in vectorized manner 
  
  k1 <- nrow(log.probs)
  k2 <- ncol(log.probs)
  u_is <- matrix(runif(k1*k2),k1,k2)
  g_is <- -log(-log(u_is))
  s <- max.col(g_is + log.probs)
  return(s) # returns a vector
}


gen.eettas <- function(di, mu, lambda, ti, kt) {
  # Generates a random number from the density p(eetta_i|the rest) in the case where 
  # the corresponding data point comes from the simulation model i.e. z_i2 == 1.
  # NOTE: It is assumed that di > 0, if di == 0 then sampling can be done straightforward way.
  # Uses the fact that the conditional density of eetta_i can be presented as finite Gamma mixture and draws exact samples. 
  # Very fast (unless d_i is large) and does not require evaluation of special functions etc. 
  # This version is a (partly) vectrorized version i.e. di and ti can be vectors.
  
  n <- length(ti)
  if (n < 1) {
    return(NULL)
  }
  beta <- 1 + lambda/(2*mu)
  log.betamuti <- log(beta*mu*ti)
  k <- rep(NA,n)
  for (i in 1:n) {
    #log.w.tilde <- cumsum(c(0,log(di[i]-(0:(di[i]-1))) - log.betamuti[i]))
    # if (kt == 1) {
    #   log.w.tilde <- cumsum(c(0, log(di[i]:1) - log.betamuti[i]))
    # } else {
    #   # this actually also handles kt == 1 case...
    #   log.w.tilde <- cumsum(c(0, log((kt + 0:(di[i]-1))*(di[i]:1)/(1:di[i])) - log.betamuti[i]))
    # }
    log.w.tilde <- cumsum(c(0, log((kt + 0:(di[i]-1))*(di[i]:1)/(1:di[i])) - log.betamuti[i]))
    
    # naive:
    #w <- exp(log.w.tilde)/sum(exp(log.w.tilde))
    #k[i] <- sample(x = 1:di[i], size = 1, replace = TRUE, prob = w)
    
    # robust:
    k[i] <- sample.cat.density(log.probs = log.w.tilde)
  }
  return(rgamma(n = n, shape = k + kt - 1, rate = 2*beta))
}


########################################################################################

gen.neff.mu.proposal <- function(neff.mu.cur, neff.v, mu.v) {
  # Generates one sample point (n_eff,mu) from a proposal density q. This is needed in the Metropolis
  # step of the Metropolis-within-Gibbs sampling algorithm. Note that n_eff is discrete while mu is 
  # continuous so this density is mixed random vector in ZxR
  # TODO: correlation between neff and mu could be taken into account
  
  if (length(neff.mu.cur) != 2) {
    stop('Incorrect neff,mu current point.')
  }
  
  # this is theoretically valid approach and this proposal is symmetric -> its pdf cancels out 
  # in Metropolis acceptance probability computation.
  neff.s <- ceiling(sqrt(neff.v))
  neff.vals <- seq(neff.mu.cur[1] - 4*neff.s, neff.mu.cur[1] + 4*neff.s, by=1) 
  neff.gaus.w <- exp(-1/(2*neff.v)*(neff.mu.cur[1] - neff.vals)^2)
  neff.gaus.w <- neff.gaus.w / sum(neff.gaus.w)
  neff.new <- sample(x = neff.vals, size = 1, prob = neff.gaus.w)
  mu.new <- rnorm(n = 1, mean = neff.mu.cur[2], sd = sqrt(mu.v))

  
  # TESTING: ALLOW OCCASIONAL LARGER JUMPS TO BE ABLE TO TRAVEL BETWEEN POSSIBLE MODES BETTER
  if (T) {
    r <- runif(1)
    if (r > 0.9) {
      nl <- 2000
      neff.new <- sample(x = seq(neff.mu.cur[1]-nl, neff.mu.cur[1]+nl, by=1), size = 1)
    } else if (r > 0.8) {
      ml <- 0.001
      mu.new <- neff.mu.cur[2] + (-ml+2*ml*runif(1))
    }
  }
  return(c(neff.new, mu.new))
}


########################################################################################

get.2x2.grid <- function(x.grid, y.grid, xy) {
  # Returns the minimal x and y index of the 2x2 grid containing the point xy.
  # This is for bilinear interpolation.
  
  # check input
  #if (xy[1] < x.grid[1] || xy[1] > x.grid[length(x.grid)] || xy[2] < y.grid[1] || xy[2] > y.grid[length(y.grid)]
  #    || length(xy) != 2) {
  #  stop('Incorrect input to determining the 2x2 grid.')
  #}
  
  # find the 2x2 area containing the point 'xy'
  x1 <- min(c(which(x.grid>xy[1]), length(x.grid))) - 1
  y1 <- min(c(which(y.grid>xy[2]), length(y.grid))) - 1
  return(c(x1,y1)) 
}


interp.bilinear <- function(f.grid, x.grid, y.grid, x1, y1, xy) {
  # Bilinear interpolation to approximate p_S(d_i1|n_eff,mu) and other functions at such n_eff,mu-points 
  # that are not in the simulation grid. 
  # NOTE: Produces continuous but not differentiable (at the edges) interpolant.
  # TODO: Use bicubic intepolation instead?
  
  # check input
  #if (xy[1] < x.grid[1] || xy[1] > x.grid[length(x.grid)] || xy[2] < y.grid[1] || xy[2] > y.grid[length(y.grid)]
  #   || length(xy) != 2 || dim(f.grid)[1] != length(x.grid) || dim(f.grid)[2] != length(y.grid)) {
  #  stop('Incorrect input to bilinear interpolation.')
  #}
  
  # approximate with the value in one of the corner points 
  #return(f.grid[x1+1,y1+1])
  
  x2 <- x1 + 1
  y2 <- y1 + 1
  fQ <- f.grid[x1:x2, y1:y2]
  h1h2 <- (x.grid[x2]-x.grid[x1])*(y.grid[y2]-y.grid[y1])
  
  # linear interpolation wrt mu
  #h2 <- y.grid[y2] - y.grid[y1]
  #return(fQ[1,1] + (xy[2] - y.grid[y1])*(fQ[1,2] - fQ[1,1])/h2)
  
  # just mean of the neighbours
  #return(mean(fQ))
  
  # bilinear interpolation
  fxy <- c(x.grid[x2]-xy[1], xy[1]-x.grid[x1]) %*% fQ %*% (c(y.grid[y2]-xy[2], xy[2]-y.grid[y1]))
  # NaN output can occur if f is -Inf at some point of the 2x2 area but we evaluate at some corner point
  fxy[is.na(fxy)] <- -Inf
  return(as.numeric(fxy/h1h2))
}


interp.bilinear.vectorized <- function(f.grid, x.grid, y.grid, x1, y1, xy) {
  # Bilinear interpolation to approximate p_S(d_i1|n_eff,mu) and other functions at such n_eff,mu-points 
  # that are not in the simulation grid. 
  # This version is vectorized wrt f.grid -- useful when need to compute many interpolations
  # at the same point xy and the same underlying grid but for different data of f values.
  # NOTE: Produces continuous but not differentiable (at the edges) interpolant.
  # TODO: Use bicubic intepolation instead?
  
  x2 <- x1 + 1
  y2 <- y1 + 1
  h1h2 <- (x.grid[x2]-x.grid[x1])*(y.grid[y2]-y.grid[y1])
  
  # bilinear interpolation
  x2x <- x.grid[x2] - xy[1]
  xx1 <- xy[1] - x.grid[x1]
  y2y <- y.grid[y2] - xy[2]
  yy1 <- xy[2] - y.grid[y1]
  if (is.matrix(f.grid)) {
    # ensures that the following works for 2d inputs as well
    dim(f.grid)[3] <- 1 
  }
  
  fxy <- f.grid[x1,y1,]*x2x*y2y + f.grid[x2,y1,]*xx1*y2y + f.grid[x1,y2,]*x2x*yy1 + f.grid[x2,y2,]*xx1*yy1
  # NaN output can occur if f is -Inf at some point of the 2x2 area but we evaluate at some corner point
  fxy[is.na(fxy)] <- -Inf
  return(fxy/h1h2)
}


test.interp <- function() {
  # Simple test case for the interpolation code.
  
  graphics.off()
  x_grid <- c(20,2000)
  y_grid <- 0:1
  data <- matrix(c(1,1,2,4),2,2)
  xy <- c(200,0)
  print(data)
  
  s <- get.2x2.grid(x_grid, y_grid, xy) 
  print(s)
  
  h1h2 <- (x_grid[2]-x_grid[1])*(y_grid[2]-y_grid[1])
  i <- interp.bilinear(data, x_grid, y_grid, s[1], s[2], xy, h1h2)
  print(i)
}


########################################################################################

get.z.plot.inds <- function(d_is) {
  # Returns example indexes of z for plotting.
  
  example_ds <- c(0,5,10,25,50,100,500,1000)
  example_ds <- example_ds[example_ds <= max(d_is)]
  inds <- NULL
  for (i in example_ds) {
    inds <- c(inds, min(which(i <= d_is)))
  }
  return(unique(inds[is.finite(inds)]))
}


simple.sample.plot <- function(zs, omegas, lambdas, neffs, mus, t_0s, unnorm.post, t_is, d_is, 
                               neff.grid, mu.grid, z.plot.inds, fig.name) {
  # Simple code for plotting (some of) the MCMC chains
  # TODO: if MCMC sampling fails, some errors can occur here
  
  library(latex2exp)
  
  max.nr.plots <- 5
  max.nr.samples <- 2000
  inds <- unique(floor(seq(1,length(mus),len=max.nr.samples)))
  inds.x <- 1:length(inds)
  
  # which of the N z and t_0 parameters are plotted?
  if (is.null(z.plot.inds)) {
    z.plot.inds <- unique(floor(seq(1,dim(zs)[1],len=max.nr.plots)))
  }
  
  hei <- 3 + length(z.plot.inds)
  pdf(file=fig.name, width=15, height=1.5*hei, pointsize=12) 
  par(mar=2*c(1,1,1,1))
  par(mfrow=c(hei,4))
  
  breaks01 <- seq(0,1,len=50)
  breaks01z <- seq(0,1,len=10)
  max.nr.breaks <- 1000
  #breaks.neff <- seq(min(neff.grid),max(neff.grid),len=min(max.nr.breaks,length(neff.grid)+1))
  #breaks.mu <- seq(min(mu.grid),max(mu.grid),len=min(max.nr.breaks,length(mu.grid)+1))
  #breaks.neff <- seq(min(neffs),max(neffs),len=min(max.nr.breaks,length(neff.grid)+1))
  #breaks.mu <- seq(min(mus),max(mus),len=min(max.nr.breaks,length(mu.grid)+1))
  breaks.neff <- seq(min(neffs),max(neffs),len=15)
  breaks.mu <- seq(min(mus),max(mus),len=15)
  
  # first row
  #plot(inds.x, neffs[inds], main = 'n_eff', ylim = c(min(neff.grid),max(neff.grid)))
  plot(inds.x, neffs[inds], main = TeX('$n_{eff}$'), ylim = c(min(neffs),max(neffs)))
  hist(neffs[inds], main = TeX('$n_{eff}$'), breaks = breaks.neff)
  #plot(inds.x, mus[inds], main = 'mu', ylim = c(min(mu.grid),max(mu.grid)))
  plot(inds.x, mus[inds], main = TeX('$\\mu$'), ylim = c(min(mus),max(mus)))
  hist(mus[inds], main = TeX('$\\mu$'), breaks = breaks.mu)
  
  # second row
  plot(inds.x, omegas[1,inds], main = TeX('$\\omega_{1}$'), ylim = c(-0.1,1.1))
  hist(omegas[1,inds], main = TeX('$\\omega_{1}$'), xlim = c(0,1), breaks = breaks01)
  plot(inds.x, omegas[2,inds], main = TeX('$\\omega_{2}$'), ylim = c(-0.1,1.1))
  hist(omegas[2,inds], main = TeX('$\\omega_{2}$'), xlim = c(0,1), breaks = breaks01)
  
  # third row
  plot(inds.x, lambdas[inds], main = TeX('$\\lambda$'))
  hist(lambdas[inds], main = TeX('$\\lambda$'))
  # computations of the log post rarely fails -> just Inf/NaN/NA
  if (any(is.finite(unnorm.post[inds]))) {
    plot(inds.x, unnorm.post[inds], main = TeX('log posterior'))
    hist(unnorm.post[inds], main = TeX('log posterior'))
  } else {
    #browser()
  }
  
  # fourth row etc., z and t_i0
  # there are N of these plots so plotting them all nicely is not possible
  for (i in z.plot.inds) {
    tdi <- paste(' (d_i=',d_is[i],', t_i=',round(t_is[i]),')',sep='') # which d_i,t_i corresponds to the plotted posterior
    
    # z-values:
    plot(inds.x, zs[i,1,inds], main = TeX(paste('$z_{', '1,', i, '}', tdi, '$', sep = '')), ylim = c(-0.1,1.1))
    barplot(table(factor(as.integer(zs[i,1,inds]),levels=c(0,1))), main = TeX(paste('$z_{', '1,', i, '}', tdi, '$', sep = '')), 
            col = 'white', axes = TRUE, space = 0)
    
    # t_0-values:
    plot(inds.x, t_0s[i,inds], main = TeX(paste('t_{0,', i, '}', tdi, '$', sep = '')))
    hist(t_0s[i,inds], main = TeX(paste('t_{0,', i, '}', tdi, '$', sep = '')), breaks = 50)
  }
  
  # some debugging if grid evaluations are used...
  nr.unique.neff <- length(unique(neffs))
  nr.unique.mu <- length(unique(mus))
  #print(paste('Number of unique sampled n_eff grid points = ', nr.unique.neff, sep=''))
  #print(paste('Number of unique sampled mu grid points = ', nr.unique.mu, sep=''))
  if (nr.unique.neff == 1) {
    warning('Only one n_eff grid value in the samples.', immediate. = TRUE)
  }
  if (nr.unique.mu == 1) {
    warning('Only one mu grid value in the samples.', immediate. = TRUE)
  }
  
  # save the figure 
  dev.off()
}


plot.neff.mu.joint.post <- function(neffs, mus, neff.grid, mu.grid, plot.grid = TRUE, fig.name = NULL, 
                                    plot.density.kde = TRUE, plot.first.point = TRUE, true.params = NULL,
                                    minimal.grid = FALSE, titles = NULL) {
  # Plots the n_eff, mu joint 2d posterior i.e. plots the 2d samples
  # neffs and mus can be now lists that contain the samples, in that cases the plots are included to one figure.
  # Grids are then assumed to be the same. 
  
  library(latex2exp)
  max.nr.samples <- 2000
  
  # This now plots the figures side by side if neffs and mus are lists that contain the samples
  if (is.list(neffs)) {
    neffs.all <- neffs
    mus.all <- mus
  } else {
    neffs.all <- list(neffs)
    mus.all <- list(mus)
  }
  n.plots <- length(neffs.all)
  
  # figure out the grid to use in the plottings
  if (minimal.grid) {
    # Find the smallest grid that contains the posterior(s) and plot this area only
    n <- c(min(unlist(neffs.all),na.rm = T), max(unlist(neffs.all),na.rm = T))
    m <- c(min(unlist(mus.all),na.rm = T), max(unlist(mus.all),na.rm = T))
    ln <- max(which(neff.grid<=n[1])); un <- min(which(neff.grid>=n[2]))
    lm <- max(which(mu.grid<=m[1])); um <- min(which(mu.grid>=m[2]))
    neff.grid <- neff.grid[ln:un]
    mu.grid <- mu.grid[lm:um]
  } # otherwise use the provided full grid
  
  if (length(fig.name) >= 1) {
    if (n.plots > 1) {
      le <- 4
      pdf(file=fig.name, width = n.plots * le, height = le, pointsize = 14)
      par(mfrow = c(1, n.plots))
      par(mai = c(0.6, 0.58, 0.3, 0.2))
    } else {
      le <- 6
      pdf(file=fig.name, width = le, height = le)
    }
  }
  
  # Draw the content in the figure
  for (i in 1:n.plots) {
    neffs.i <- neffs.all[[i]]
    mus.i <- mus.all[[i]]
    inds <- unique(floor(seq(1,length(mus.i),len=max.nr.samples)))
    
    # plot simulated points
    plot(neffs.i[inds], mus.i[inds], type = 'p', pch=19, cex=0.25, col = 'blue', xaxs="i", yaxs="i",
         xlab = TeX('$n_{eff}$'), ylab = TeX('$\\mu$'), cex.lab = 1.3,
         xlim = c(min(neff.grid),max(neff.grid)), ylim = c(min(mu.grid),max(mu.grid)))
    if (!is.null(titles)) {
      title(main = TeX(titles[[i]])) # plot title if it is given as input
    }
    grid()
    
    # plot the 2d density estimate by kde
    if (plot.density.kde) {
      library('MASS')
      nl <- length(neff.grid)
      ml <- length(mu.grid)
      max.nr.samples <- 2000
      inds <- unique(floor(seq(1,length(mus.i),len=max.nr.samples)))
      f1 <- kde2d(neffs.i[inds], mus.i[inds], n = c(4*nl,4*ml), 
                  lims = c(min(neff.grid),max(neff.grid),min(mu.grid),max(mu.grid)))
      contour(f1, drawlabels = FALSE, nlevels = 15, add = TRUE)
    }
    
    # plot the initial point of the chain
    if (plot.first.point) {
      lines(neffs.i[1], mus.i[1], pch = 24, cex = 1.4, col = 'red', bg = 'red', type = 'p') 
    }
    
    # if artificial data was used and true param value was given as input, plot it also
    if (!is.null(true.params) && !is.null(true.params[['neff']]) && !is.null(true.params[['mu']])) {
      lines(true.params$neff, true.params$mu, pch = 23, cex = 1.4, col = 'green', bg = 'green', type = 'p')
    }
    
    # plot grid points
    if (plot.grid) {
      xy <- meshgrid(neff.grid, mu.grid)
      points(as.vector(xy$x), as.vector(xy$y), pch = '.', col='grey')
    }
  }
  
  # save the figure if filename was given as input
  if (length(fig.name) >= 1) {
    dev.off()
  }
}




