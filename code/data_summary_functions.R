# Functions for computing distances, summaries and the discrepancy needed for the ABC inference.

count.total.number.mutations <- function(pop) {
	
	pop$snp.sets
	
	snp.set.sizes <- unlist(lapply(pop$snp.sets, length))
	num.snps.in.strains <- snp.set.sizes[pop$strains]
	total.number.mutations <- sum(num.snps.in.strains)
	
	return(total.number.mutations)
	
}


compute.distances.in.simulation <- function(save.path, save.extension, generations.to.consider, dist.path) {
	
	for (generation.index in generations.to.consider) {
		
		population.file.name <- paste(save.path, '/pop', generation.index, save.extension, '.RData', sep='')
		
		load(population.file.name) # pop
		
		distance.distribution <- compute.distance.distribution(pop=pop)
		
		distance.file.name <- paste(dist.path, '/distances', generation.index, save.extension, '.RData', sep='')
		
		save(distance.distribution, file=distance.file.name) 
		
	}
	
}


compute.distance.distribution <- function(pop) {
	# The value returned by this function is a list with two
	# fields:
	# 
	# "snp.set.dist" is a distance matrix for SNP sets
	# (the number of SNPs by which two sets differ)
	# 
	# "dist.elem.count" is also a matrix. Element (i,j)
	# tells how many strain pairs there are, such that the
	# first strain has SNP set i, and the second has set j.
	#
	# So, together these two specify the full distance
	# distribution between all strain pairs. The distance
	# between two specific strains x and y can be obtained 
	# as snp.set.dist[pop$strains[x],pop$strains[y]]
	
	# Compute distance matrix between snp.sets
	n.snp.sets <- length(pop$snp.sets)
	snp.set.dist <- matrix(0, n.snp.sets, n.snp.sets)
	dist.elem.count <- matrix(0, n.snp.sets, n.snp.sets)
	
	strain.type.counts <- table(pop$strains)
	n.within.pairs <- choose(strain.type.counts,2)
	# Check that all SNP sets are present in the population
	if (!all(as.numeric(names(strain.type.counts)) == 1:n.snp.sets)) {
		stop('Corrupt data structure')
	}
	
	# Note the special case of no mutations: then snp.sets contain only empty list, the code now handles this
	if (n.snp.sets == 1) {
	  snp.set.dist <- matrix(0, 1, 1)
	  dist.elem.count <- matrix(0, 1, 1)
	  dist.elem.count[1, 1] <- n.within.pairs["1"]
	}
	else if (n.snp.sets > 1) {
	  for (set1.index in seq(1, n.snp.sets - 1)) {
	    set1 <- pop$snp.sets[[set1.index]]
	    
	    for (set2.index in seq(set1.index + 1, n.snp.sets)) {
	      set2 <- pop$snp.sets[[set2.index]]
	      
	      # Distance is the number of SNPs that are found
	      # in exactly one or the other of the SNP sets,
	      # but not in both.
	      snp.set.dist[set1.index, set2.index] <- length(union(set1, set2)) - length(intersect(set1, set2))
	      
	      # Compute the number of strain pairs that
	      # have this distance, i.e., one of the strains
	      # has set1 and the other has set2.
	      dist.elem.count[set1.index, set2.index] <- strain.type.counts[as.character(set1.index)] * strain.type.counts[as.character(set2.index)]
	      
	    }
	    dist.elem.count[set1.index, set1.index] <- n.within.pairs[as.character(set1.index)]
	  }
	  dist.elem.count[n.snp.sets, n.snp.sets] <- n.within.pairs[as.character(n.snp.sets)]
	} else {
	  stop('Corrupt data structure')
	}
	
	res <- list()
	res$snp.set.dist <- snp.set.dist + t(snp.set.dist) # Symmetric
	res$dist.elem.count <- dist.elem.count
	res$gen <- pop$gen
	
	return(res)
}


############################################
# Summaries and discrepancy for ABC fitting
############################################

summaries.one.patient <- function(distance.distribution) {
  # Computes the summary vector for one distance distribution (one patient)
  
  # compute the average distance over all the pairwise distances at that time point
  n.pairs <- sum(distance.distribution$dist.elem.count)
  d.sum <- distance.distribution$snp.set.dist * distance.distribution$dist.elem.count
  s1 <- sum(d.sum)/n.pairs
  
  # if other summaries are also used, add them here

  return(s1)
}


discrepancy <- function(discrepancy.type, dists, dists.AM, data, data.AM, save.interval, save.interval.AM, 
                        print.discrepancy=FALSE, neff=NULL, mu=NULL) {
  # Computes the discrepancy for fitting the simulation model with ABC, all data
  # points (i.e. patients 1-8 and A-M) can be used for computing the discrepancy. 
  
  rem.index <- 2 # index of patient #1219 in the data
  
  # NOTE: THIS NOW WORKS ONLY IF LENGTH(SAVE.INTERVAL) > 1
  if (length(save.interval) <= 1) {
    stop('This now works only if save time points are given i.e. save.interval must be a vector longer than 1.')
  }
  
  # data lengths
  n.first <- length(data$first)
  n.AM <- length(data.AM)
  n.s.AM <- length(save.interval.AM)
  
  ## DISCREPANCY "1"
  if (discrepancy.type$type == 'l1_sep') {
    ## Compute the L1 distances separately and then the total discrepancy
    # patients 1-8
    d.first <- rep(0,n.first)
    d.last <- rep(0,n.first)
    for (i in 1:n.first) {
      # get the time of making the observation
      time1.ind <- which(discrepancy.type$times18_1[i]==save.interval)
      time2.ind <- which(discrepancy.type$times18_2[i]==save.interval)
      if (length(time1.ind) < 1 || length(time2.ind) < 1) {
        stop('Time of observed data was not simulated.')
      }
      d.first[i] <- l1.distance.discrete(dists[[i]][[time1.ind]], data$first[[i]])
      d.last[i] <- l1.distance.discrete(dists[[i]][[time2.ind]], data$last[[i]])
    }
    
    # deal with patient #1219
    if (discrepancy.type$ignore1219) {
      d.first <- d.first[-rem.index]
      d.last <- d.last[-rem.index]
    }
    
    # patients A-M
    d.AM <- rep(0,n.AM)
    if (!discrepancy.type$ignoreAM) {
      for (i in 1:n.AM) {
        dis <- rep(0,n.s.AM)
        for (j in 1:n.s.AM) {
          dis[j] <- l1.distance.discrete(dists.AM[[i]][[j]], data.AM[[i]])
        }
        # Since the obs.time is not known, we take such simulated times that yield minimum discrepancy
        d.AM[i] <- min(dis)
      }
    }
    
    # total discrepancy
    d <- (sum(d.first) + sum(d.last) + sum(d.AM))/(2*length(d.first) + (!discrepancy.type$ignoreAM)*n.AM)
    if (print.discrepancy) {
      print.discrepancy.debug.l1(neff, mu, d, d.first, d.last, d.AM) 
    }
    return(d)
  }
  
  
  ## ALTERNATIVE FUNCTIONS FOR THE DISCREPANCY ARE ADDITIONALLY COMPUTED BELOW
  
  ## 0) Compute true data summaries (this and some other values could be precomputed but it is very fast anyway)
  s.obs1 <- sapply(data$first, mean)
  s.obs2 <- sapply(data$last, mean)
  s.obs.AM1 <- sapply(data.AM, mean)

  ## Compute the summaries from the simulation model 
  # 1) Summaries for the first data set (patients 1-8, possibly #1219 to be excluded)
  s.model1 <- rep(NA,n.first)
  s.model2 <- rep(NA,n.first)
  for (i in 1:n.first) {
    time1.ind <- which(discrepancy.type$times18_1[i]==save.interval)
    time2.ind <- which(discrepancy.type$times18_2[i]==save.interval)
    if (length(time1.ind) < 1 || length(time2.ind) < 1) {
      stop('Time of observed data was not simulated.')
    }
    s.model1[i] <- summaries.one.patient(dists[[i]][[time1.ind]])
    s.model2[i] <- summaries.one.patient(dists[[i]][[time2.ind]])
  }
  
  if (discrepancy.type$ignore1219) {
    # remove the index 2 that corresponds to the patient #1219
    rem.index <- 2 
    s.obs1 <- s.obs1[-rem.index]
    s.obs2 <- s.obs2[-rem.index]
    s.model1 <- s.model1[-rem.index]
    s.model2 <- s.model2[-rem.index]
  }
  
  # 2) Summaries for the second data set (patients A-M)
  s.model.AM1 <- matrix(NA,n.s.AM,n.AM)
  for (i in 1:n.AM) {
    for (j in 1:n.s.AM) {
      s.model.AM1[j,i] <- summaries.one.patient(dists.AM[[i]][[j]])
    }
  }
  
  ## Compute the discrepancy using the summaries computed above
  d.AM <- 0
  if (discrepancy.type$type == 'eucl') {
    d1 <- norm_vec(s.obs1 - s.model1)
    d2 <- norm_vec(s.obs2 - s.model2)
    #d2 <- norm_vec((s.obs1 - s.obs2) - (s.model1 - s.model2))
    
    if (!discrepancy.type$ignoreAM) {
      abs.errors.AM <- apply(abs(t(t(s.model.AM1) - s.obs.AM1)), 2, min)
      #print('A-M abs. errors:'); print(abs.errors.AM); print(' ')
      d.AM <- norm_vec(abs.errors.AM) 
    }
    
    ## Final discrepancy over all data
    #d <- sqrt(1*d1^2 + 1*d2^2 + 1*d.AM^2)
    d <- sum(d1 + d2 + d.AM)
    
  } else {
    stop('Incorrect discrepancy type.')
  }
  
  ## Debug printing
  if (print.discrepancy) {
    print.discrepancy.debug(neff, mu, d, d1, d2, d.AM, s.obs1, s.obs2, s.obs.AM1, s.model1, s.model2, s.model.AM1)
  }
  return(d)
}


l1.distance.discrete <- function(distance.distribution1, obs.samples) {
  # Compute the l1 distance between two pmf's that are represented by samples, 
  # the domain is assumed to be non-negative integers
  
  maxd <- max(distance.distribution1$snp.set.dist, obs.samples)
  counts1 <- rep(0,maxd+1)
  counts2 <- rep(0,maxd+1)
  for (i in 0:maxd) {
    ind1 <- which(distance.distribution1$snp.set.dist == i)
    if (length(ind1) > 0) {
      counts1[i+1] <- sum(distance.distribution1$dist.elem.count[ind1])
    }
    counts2[i+1] <- sum(obs.samples == i)
  }
  return(0.5*sum(abs(counts1/sum(counts1) - counts2/sum(counts2))))
}


norm_vec <- function(x) sqrt(sum(x^2))
#norm_vec <- function(x) norm_abs_sum(x)

norm_abs_sum <- function(x) sum(abs(x))


print.discrepancy.debug.l1 <- function(neff, mu, d, d.first, d.last, d.AM) {
  # Prints some information about the discrepancy. Intended for testing. 
  
  # print parameter and discrepancy values
  cat(paste('neff: ', neff, ', mu: ', mu, ', d = ', d, sep = ''))
  cat('\n')
  
  # print discrepancy contributions from data sets
  print('data_1-8_first:')
  print(d.first)
  cat('\n')
  print('data_1-8_second:')
  print(d.last)
  cat('\n')
  print('data_AM:')
  print(d.AM)
  cat('\n\n')
}


print.discrepancy.debug <- function(neff, mu, d, d1, d2, d.AM, s.obs1, s.obs2, s.obs.AM1, s.model1, s.model2, s.model.AM1) {
  # Prints some information about the summaries and the discrepancy. Intended for testing. 
  
  # print parameter and discrepancy values
  cat(paste('neff: ', neff, ', mu: ', mu, ', d = ', d, sep = ''))
  cat('\n')
  cat(paste('d1: ', d1, ', d2: ', d2, ', d.AM: ', d.AM, sep = ''))
  cat('\n')
  
  # print summaries
  s1 <- cbind(s.obs1, s.model1)
  s2 <- cbind(s.obs2, s.model2)
  #s.AM <- cbind(s.obs.AM1, s.model.AM1)
  
  print('data_1-8_first:')
  print(s1)
  cat('\n')
  print('data_1-8_second:')
  print(s2)
  cat('\n')
  print('data_AM:')
  print(s.obs.AM1)
  print(s.model.AM1)
  cat('\n\n')
}




