### Modified function to run only kinships within each pop AND within each sp within each pop
### i.e., accounts for different species that have the exact same site name

individual_kinship_by_pop_sp <- function(dart_data, basedir, species, dataset, pop, sp, maf=0.05, mis=0.2, as_bigmat=TRUE) {
  require(SNPRelate)

  # combine species + population as group ID
  spop <- paste(sp, pop, sep="__")
  spop_levels <- unique(spop)
  
  kinlist <- list()
  nsamp <- nrow(dart_data$gt)
  igds_file <- dart2gds(dart_data, basedir, species, dataset)
  igds <- snpgdsOpen(igds_file)
  
  for (grp in spop_levels) {
    isamps <- which(spop == grp)
    
    iout <- snpgdsIBDMoM(
      igds,
      sample.id=rownames(dart_data$gt)[isamps],
      maf=maf,
      missing.rate=mis,
      num.thread=1,
      kinship=TRUE
    )
    
    ikout <- iout$kinship
    rownames(ikout) <- rownames(dart_data$gt)[isamps]
    colnames(ikout) <- rownames(dart_data$gt)[isamps]
    kinlist[[ grp ]] <- ikout
  }
  
  snpgdsClose(igds)
  
  if (as_bigmat) {
    bigmat <- matrix(0, nsamp, nsamp)
    rownames(bigmat) <- ""
    colnames(bigmat) <- ""
    
    istart <- 1
    for (i in seq_along(kinlist)) {
      im <- kinlist[[i]]
      iN <- nrow(im)
      
      if (istart == 1) {
        istop <- iN
      } else {
        istop <- istop + iN
      }
      
      bigmat[istart:istop, istart:istop] <- im
      rownames(bigmat)[istart:istop] <- rownames(im)
      colnames(bigmat)[istart:istop] <- rownames(im)
      istart <- istart + iN
    }
    
    return(bigmat)
  }
  
  return(kinlist)
}

