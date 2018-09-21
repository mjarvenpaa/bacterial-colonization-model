# Misc. functions

formulate.save.extension <- function(variables=c('mut'=0.0019, 'num'=500)) {
	
	save.extension <- ''
	var.names <- names(variables)
	for (i in 1:length(variables)) {
		name.now <- var.names[i]
		value.now <- signif(variables[i],digits=6)
		save.extension <- paste(save.extension, '_', name.now, value.now, sep='')
	}
	
	return(save.extension)
}


get.current.ids <- function(run.id, max.run.id, tot.sims) {
  # Returns the indexes corresponding to the current id
  # This is needed for splitting the computations to smaller tasks for computer cluster.
  #
  # run.id - current index for parallel runs
  # max.run.id - maximum index of the parallel runs
  # tot.sims - total number of simulations to be performed (on all nodes run in parallel)
  
  # some input checks
  if (run.id <= 0 || max.run.id <= 0 || tot.sims <= 0 || max.run.id > tot.sims) {
    stop('Incorrect input for selecting indexes.')
  }
  if (run.id > max.run.id) {
    return(NULL)
  }
  
  # Current simulation indexes corresponding to given run.id
  if (tot.sims %% max.run.id != 0) {
    stop('Current implementation requires that each parallel run contains the same amount of evaluations.')
  } else {
    le <- ceiling(tot.sims / max.run.id)
  }
  
  cur.sim.ids <- ((run.id - 1)*le + 1):(run.id*le)
  return(cur.sim.ids)
}


get.nr.individuals.from.nr.pairs <- function(nr.pairs, print.output=FALSE) {
  # Solves for nr.individuals from the equation choose(nr.individuals, 2) == nr.pairs
  # (Could be more efficiently solved using bisection search but here nr.indivs is not large so this code will suffice.)
  # (In our data case, we could also just set these also manually...)
  
  max.nr.indivs <- 30
  
  max.nr.pairs <- choose(max.nr.indivs,2)
  tries <- 1:max.nr.indivs
  
  nr.individuals <- rep(NA,length(nr.pairs))
  for (i in 1:length(nr.pairs)) {
    indi <- which(choose(tries,2)==nr.pairs[i])
    if (length(indi) == 0) {
      # either impossible case or nr.pairs too high for this code (which should not happen)
      if (nr.pairs[i] > max.nr.pairs) {
        warning('Too many pairs, this code does not handle that.')
      } else {
        warning('Impossible amount of pairs.')
      }
    } else {
      nr.individuals[i] <- tries[indi]
    }
  }
  if (print.output) {
    print(nr.individuals)
  }
  return(nr.individuals)
}


get.int.grid <- function(start, end, n, log.scale=FALSE) {
  # Get an integer interval, the usual or log-scale
  
  if (start > end) {
    stop('Starting value greater than end value.')
  } 
  if(log.scale) {
    values <- pmin(ceiling(exp(seq(log(start), log(end), len=n))), end) # unif in log-scale
  } else {
    values <- ceiling(seq(start, end, len=n)) 
  }
  return(unique(values))
}


meshgrid <- function(x,y) {
  # Similar to matlab's meshgrid
  xo <- outer(x, y, function(x,y){x})
  yo <- outer(x, y, function(x,y){y})
  return(list(x=xo, y=yo))
} 


filled.contour3 <- function (x = seq(0, 1, length.out = nrow(z)),
                             y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                             ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                             levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                             col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                             key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                             axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page
  
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  # on.exit(par(par.orig))
  # w <- (3 + mar.orig[2]) * par("csi") * 2.54
  # par(las = las)
  # mar <- mar.orig
  plot.new()
  # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}




