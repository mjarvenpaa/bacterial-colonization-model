# Some misc. plotting functions for the model.
# Functions for some other plottings etc. are elsewhere.

plot.j.clin.microbiol.distributions <- function(root) {
  # Plots the real data D_0 (a histogram for each case)
  
  data.path <- paste(root, '/data', sep='')
  rdata.file.name <- paste(data.path, '/J_CLIN_MIROBIOL_fig1.RData', sep='')
  
  load(rdata.file.name)
  
  figure.path <- paste(root, '/figures', sep='')
  
  fig.file.name <- paste(figure.path, '/j_clin_microbiol_fig1.pdf', sep='')
  
  names.to.plot <- names(table(data.mat[,'subject_time']))
  n.plots <- length(names.to.plot)
  
  n.rows.in.plot <- ceiling(n.plots/2)
  
  #browser()
  #png(file=fig.file.name, width = 7, height = 14, units = "in", pointsize = 10, bg = "white", res=300)
  pdf(file=fig.file.name, width = 5, height = n.rows.in.plot*2, pointsize = 10, bg = "white")
  par(mfrow=c(n.rows.in.plot,2))
  
  for (plot.index in 1:n.plots) {
    rows.now <- which(data.mat[,'subject_time'] == names.to.plot[plot.index])
    distances.now <- as.matrix(data.mat[rows.now, 'dist'])
    distances.now <- as.vector(distances.now)
    hist(distances.now, breaks=0:28, main=names.to.plot[plot.index], xlab='Distance')
  }
  graphics.off()
}


plot.distance.evolution <- function(save.extension, generations.to.consider, dist.path, figure.path) {
  # Plots the distance evolution when the distances are saved to individual files
  
  figure.file.name <- paste(figure.path, '/dist_evolution', save.extension, '.pdf', sep='')
  n.plots <- length(generations.to.consider)
  
  pdf(file=figure.file.name, width = 5, height = 2.7, pointsize = 10, bg = "white")
  #orig.par <- par(mfrow=c(dim1,dim2), mar=c(5,5,1,2)+0.1, oma=c(0,1,3.5,3))
  #png(file=figure.file.name, width = 8, height = 4, units = "in", pointsize = 10, bg = "white", res=300)
  
  for (gen.index in generations.to.consider) {
    
    distance.file.name <- paste(dist.path, '/distances', gen.index, save.extension, '.RData', sep='')
    load(distance.file.name)
    # Variable "distance.distribution", which has two fields:
    # "snp.set.dist" and "dist.elem.count"
    plot.distance.distribution(distance.distribution, gen.index)
  }
  
  #par(orig.par)
  graphics.off()
}


plot.distance.evolution.pops <- function(dists, generations.to.consider, save.extension, figure.path, adapt.bins = FALSE) {
  # Plots the distance evolution when the distances are given as input in a list datastructure
  
  # find the indexes for generations.to.consider
  pops.gens <- unlist(lapply(dists, '[[', 'gen'))
  generations.to.consider <- unique(intersect(generations.to.consider,pops.gens))
  gen.inds <- match(generations.to.consider,pops.gens)
  
  figure.file.name <- paste(figure.path, '/dist_evolution', save.extension, '.pdf', sep='')
  n.plots <- length(gen.inds)
  
  pdf(file=figure.file.name, width = 5, height = 2.7, pointsize = 10, bg = "white")
  #orig.par <- par(mfrow=c(dim1,dim2), mar=c(5,5,1,2)+0.1, oma=c(0,1,3.5,3))
  #png(file=figure.file.name, width = 8, height = 4, units = "in", pointsize = 10, bg = "white", res=300)
  
  for (gen.index in gen.inds) {
    
    # Variable "distance.distribution", which has two fields:
    # "snp.set.dist" and "dist.elem.count"
    distance.distribution <- dists[[gen.index]]
    generation.index <- dists[[gen.index]]$gen
    plot.distance.distribution(distance.distribution, generation.index=generation.index, adapt.bins=adapt.bins)
  }
  
  #par(orig.par)
  graphics.off()
}



plot.distance.distribution <- function(distance.distribution, generation.index=0, max.dist=30, adapt.bins=FALSE) {
  # A (binned) barplot of the distances
  #
  # NOTE: generation.index is used only in the plot title
  
  bin.counts <- compute.distance.histogram(distance.distribution, max.dist, adapt.bins)
  
  barplot(bin.counts, main=paste('Generation ', generation.index, sep=''), xlab="SNP distance", space=0)
  axis(1)
}


compute.distance.histogram <- function(distance.distribution, max.dist=30, adapt.bins=FALSE) {
  # Computes the distance histogram i.e. the amount of elements in each bin
  #
  # NOTE: If adapt.bins is TRUE, then the bins are set according to the data points and can vary between different data sets 
  
  if (adapt.bins) {
    bin.lb <- 0
    #bin.lb <- floor(min(distance.distribution$snp.set.dist))
    bin.ub <- max(max.dist,ceiling(max(distance.distribution$snp.set.dist))) + 1
    bin.ub <- min(bin.ub,100)
    bin.boundaries <- seq(bin.lb, bin.ub, by=1)
  } else {
    bin.boundaries <- seq(0, max.dist, by=1)
  }
  
  n.bins <- length(bin.boundaries) - 1
  bin.counts <- rep(0, n.bins)
  
  for (bin.index in 1:n.bins) {
    lb <- bin.boundaries[bin.index]
    ub <- bin.boundaries[bin.index + 1]
    elements.of.this.length <- which((distance.distribution$snp.set.dist >= lb) & (distance.distribution$snp.set.dist < ub))
    bin.counts[bin.index] <- sum(distance.distribution$dist.elem.count[elements.of.this.length])
  }
  
  # check if there are distances that exceed the highest bin boundary
  if (any(distance.distribution$snp.set.dist >= ub)) {
    warning("Distance(s) larger than max.dist were detected.")
  }
  
  return(bin.counts)
}




