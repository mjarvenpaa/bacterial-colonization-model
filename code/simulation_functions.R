# Functions for running the W-F simulation model. The simulated populations are either saved to file immediately once computed
# or, alternatively, returned in a list datastructure to be analysed elsewhere. 

simulate.evolution <- function(mutation.rate, n.strains, genome.length, n.generations.to.simulate, save.interval, save.path, 
                               save.extension, print.output=TRUE, save.output=TRUE) {
  # Main function for running the Wright-Fisjer simulations.
  # NOTE: Modified the original code so that save.interval is interpreted as a collection of time points if its length is > 1. 
  # If the length is 1, then it is interpreted as a time interval for saving (as before).
  
  time.start <- proc.time()['elapsed']
  
  pop <- initialize.population(n.strains=n.strains, genome.length=genome.length)
  
  pops <- list()
  
  for (generation.index in 1:n.generations.to.simulate) {
    
    
    # Sample the next generation
    pop <- sample.next.generation(pop=pop)
    
    
    # Add mutations
    pop <- add.mutations(pop=pop, mutation.rate=mutation.rate)
    
    
    # Remove SNP sets which none of the strains has anymore
    # as a consequence of sampling and new mutations.
    pop <- remove.unnecessary.snp.sets(pop=pop)
    
    
    # save the information that which generation is currently
    pop$gen <- generation.index
    
    
    # Save at regular intervals (if save.interval is one number) 
    # or only at specified generations (if save.interval is vector with length > 1)
    if ((length(save.interval) == 1 && (generation.index %% save.interval == 0 
                                        || generation.index == n.generations.to.simulate))
        || (length(save.interval) > 1 && any(generation.index == save.interval))) {
      
      elapsed <- proc.time()['elapsed'] - time.start
      pop$elapsed <- elapsed
      
      if (save.output && length(save.interval) == 1) {
        # population is saved to file
        save(pop, file=paste(save.path, '/pop', generation.index, save.extension, '.RData', sep=''))
      } else {
        # pop data is saved to a datastructure that is returned instead of always saving pop data to a new file
        # (which can result huge number of small files!)
        pops[[length(pops)+1]] <- pop
      }
      
      if (!check.population.validity(pop)) {
        stop('Corrupt population structure')
      }
      
      if (print.output) {
        print(paste('Gen ', generation.index, ' complete. Time=', signif(elapsed, digits = 3), '.', sep=''))
      }
    }
  }
  return(pops)
}



initialize.population <- function(n.strains, genome.length) {
  
  # The idea in data structures is the following. pop$snp.sets is a list of 
  # sets of SNPs  observed in some genome in the population. pop$strains maps
  # each strain to one of these sets. The idea is that this way the sampling
  # of generations can be done very rapidly, by only updating elements of in
  # a vector rather than rows in a matrix.
  
  pop <- list()
  pop$snp.sets <- list()
  pop$snp.sets[[1]] <- numeric() # Genome with no SNPs.
  # Assumptions:
  # 1) SNPs in each SNP set are ordered all the time.
  # 2) Different SNP sets must be truly different.
  
  # In the beginning there are only one kind of 
  # strains, corresponding to the first SNP set.
  pop$strains <- rep(1, n.strains)
  
  pop$n.strains <- n.strains
  
  pop$genome.length <- genome.length
  
  return(pop)
}



sample.next.generation <- function(pop) {	
  
  n.strains <- pop$n.strains
  
  offspring.indices <- sample(x=n.strains, size=n.strains, replace=T)	
  
  pop$strains <- pop$strains[offspring.indices]
  
  return(pop)
}



remove.unnecessary.snp.sets <- function(pop) {
  # Remove such snp sets that none of the strains has.
  # Update the snp set indices correspondingly.
  
  max.set.index <- length(pop$snp.sets)
  snp.sets.present <- unique(pop$strains)
  
  unnecessary.sets <- setdiff((1:max.set.index), snp.sets.present)
  
  if (length(unnecessary.sets) > 0) {
    to.keep <- 1:max.set.index
    to.keep <- to.keep[-unnecessary.sets]
    mapping.old.to.new <- match(1:max.set.index, to.keep)
    # This tells what the new set index of each remaining 
    # set will be after the removal of the unnecessary sets.
    
    pop$strains <- mapping.old.to.new[pop$strains]
    
    pop$snp.sets <- pop$snp.sets[-unnecessary.sets]
  }
  
  return(pop)
}



add.mutations <- function(pop, mutation.rate) {
  # mutation.rate is the expected number of 
  # mutations per strain per generation.
  
  n.strains <- length(pop$strains)
  
  number.of.possible.sites.to.mutate <- n.strains * pop$genome.length
  
  num.mutations <- rbinom(n=1, size=number.of.possible.sites.to.mutate, prob=mutation.rate/pop$genome.length)
  
  if (num.mutations > 0) {
    
    elements.to.mutate <- sample(pop$genome.length * n.strains, num.mutations, replace=FALSE)
    # We assume that each mutation happens in a different position,
    # justified by the negligible probability of two mutations occurring
    # in the same position in the same strain at once. The previous line
    # gives the indexes of the elements in the n.strains*genome.length table
    # that will be mutated.
    
    strains.to.mutate <- elements.to.mutate %% n.strains
    zero.strains <- which(strains.to.mutate == 0)
    if (length(zero.strains) > 0) {
      strains.to.mutate[zero.strains] <- n.strains
    }
    genome.positions.to.mutate <- ceiling(elements.to.mutate / n.strains)
    
    
    for (mutation.index in 1:num.mutations) {
      # First compute the updated SNP set
      strain.index <- strains.to.mutate[mutation.index]
      position.to.mutate <- genome.positions.to.mutate[mutation.index]
      
      old.snp.set <- pop$snp.sets[[pop$strains[strain.index]]]
      if (length(old.snp.set)==0) {
        
        new.snp.set <- position.to.mutate
        
      } else if (any(old.snp.set==position.to.mutate)) {
        
        # Remove the SNP from list
        if (length(old.snp.set)==1) {new.snp.set <- numeric()} 
        else {new.snp.set <- sort(setdiff(old.snp.set, position.to.mutate))}
        
      } else {
        
        # Just add the position to the old SNP set
        new.snp.set <- sort(union(old.snp.set, position.to.mutate))
        
      }
      
      # Check if the updated SNP set is equal
      # to some already existing SNP set
      # THIS SEEMS TO TAKE MOST OF THE COMPUTATION TIME
      is.same.set <- lapply(pop$snp.sets, function(x){
        return.value <- FALSE
        if (length(x)==length(new.snp.set)) {
          if (all(x==new.snp.set)) {
            # For this to work, SNP sets must be ordered
            return.value <- TRUE
          }
        }
        return(return.value)
      })
      is.same.set <- unlist(is.same.set) # Max one should be true
      
      
      same.set.index <- which(is.same.set)
      if (length(same.set.index)==0) {
        # The new SNP set did not exist before
        
        # Select a label for the new SNP set
        new.snp.set.label <- length(pop$snp.sets) + 1
        
        # Update the strain's pointer to the new SNP set
        pop$strains[strain.index] <- new.snp.set.label
        
        # Update SNP set descriptions
        pop$snp.sets[[new.snp.set.label]] <- new.snp.set
        
      } else if (length(same.set.index)==1) {
        # The updated SNP set is one of the already existing SNP sets.
        # Just update the strain's pointer to the SNP set.
        pop$strains[strain.index] <- same.set.index
        
      } else {
        save(list=ls(), file=paste('error_state.RData', sep=''))
        stop('Multiple identical SNP sets')
      }
    }
  }
  
  return(pop)
}


check.population.validity <- function(pop) {
  
  is.valid <- TRUE
  
  # Check that SNP sets are ordered
  n.snp.sets <- length(pop$snp.sets)
  for (i in 1:n.snp.sets) {
    snp.set.now <- pop$snp.sets[[i]]
    if (length(snp.set.now) > 1) {
      sorted.set <- sort(snp.set.now)
      if (!all(sorted.set == snp.set.now)) {
        is.valid <- FALSE
      }
    }
  }
  #print(is.valid)
  
  # Check that SNP sets are distinct, but only if there are more than 1 set.
  if (n.snp.sets > 1) {
    for (set1 in seq(1, n.snp.sets-1)) {
      for (set2 in seq(set1+1, n.snp.sets) ) {
        snps1 <- pop$snp.sets[[set1]]
        snps2 <- pop$snp.sets[[set2]]
        #tryCatch(snps2 <- pop$snp.sets[[set2]], error = browser(), finally = {})
        if (length(snps1) == length(snps2)) {
          
          if (length(snps1) == 0) {
            # There are two sets with zero elements
            is.valid <- FALSE
          }
          
          if (all(snps1 == snps2)) {
            is.valid <- FALSE
          }
        }
      }
    }
  }
  #print(is.valid)
  
  # Check that there are no extra SNP-sets that no strain has
  snp.sets.present <- sort(unique(pop$strains))
  if (length(snp.sets.present) != n.snp.sets) {
    is.valid <- FALSE
  } else if (!all(snp.sets.present == seq(1,n.snp.sets))) {
    is.valid <- FALSE
  }
  #print(is.valid)
  
  return(is.valid)
}


resample.pop <- function(pops, n.strains.new) {
  # Select a smaller subset of strains and update the populations (pop) accordingly
  #
  # This is needed for computing the discrepancy since the true observed data has been measured
  # from a small subset of the genomes and not using all of them. Thus there is additional
  # stochasticity due to the variation of selecting the genomes to be measured.
  
  for (i in 1:length(pops)) {
    
    n.strains <- pops[[i]]$n.strains
    if (n.strains.new > n.strains) {
      stop('Error in population resampling step. Observed population n.strains larger than simulated one.')
    }
    
    # resample genomes
    offspring.indices <- sample(x=n.strains, size=n.strains.new, replace=F)	
    
    pops[[i]]$strains <- pops[[i]]$strains[offspring.indices]
    
    # update also the size value in each pop
    pops[[i]]$n.strains <- n.strains.new
    
    # Remove SNP sets which none of the strains has anymore
    # as a consequence of sampling.
    pops[[i]] <- remove.unnecessary.snp.sets(pop=pops[[i]])
    
    if (!check.population.validity(pops[[i]])) {
      stop('Corrupt population structure in resampling step.')
    }
  }
  
  return(pops)
}


random.pop.distance <- function(pop, n = 1) {
  # Selects two strains randomly from population 'pop' and returns their distance.
  # This is repeated 'n' times.
  
  distances <- rep(NA,n)
  n.strains <- pop$n.strains
  
  for (i in 1:n) {
    # select 2 strains
    two.strains <- sample(x=n.strains, size=2, replace=FALSE)
    
    # compute their distance, as in an another function elsewhere 
    # except that here we do not need to compute the whole distance matrix
    set1 <- pop$snp.sets[[pop$strains[two.strains[1]]]]
    set2 <- pop$snp.sets[[pop$strains[two.strains[2]]]]
    distances[i] <- length(union(set1, set2)) - length(intersect(set1, set2))
  }
  return(distances)
}


generate.toy.population.for.testing <- function() {
  # This is only for testing
  
  pop <- list()
  
  pop$strains <- c(5,2,1,1,4)
  
  pop$snp.sets <- list()
  pop$snp.sets[[1]] <- c(5,10,12)
  pop$snp.sets[[2]] <- numeric()
  pop$snp.sets[[3]] <- c(1)
  pop$snp.sets[[4]] <- c(1,5)
  pop$snp.sets[[5]] <- c(3,5)
  pop$snp.sets[[6]] <- c(1,5,10,12)
  
  return(pop)
}





