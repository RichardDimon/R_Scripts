### Modified function to run only kinships within each pop AND within each sp within each pop
### i.e., accounts for different species that have the exact same site name

individual_kinship_by_pop_sp <- function(dart_data, basedir, species, dataset, pop, sp, maf=0.1, mis=0.2, as_bigmat=TRUE) {
  require(SNPRelate)
  
  # create species-pop ID
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
      sample.id = rownames(dart_data$gt)[isamps],
      maf = maf,
      missing.rate = mis,
      num.thread = 1,
      kinship = TRUE
    )
    
    ikout <- iout$kinship
    rownames(ikout) <- rownames(dart_data$gt)[isamps]
    colnames(ikout) <- rownames(dart_data$gt)[isamps]
    
    kinlist[[grp]] <- ikout
  }
  
  snpgdsClose(igds)
  
  if (!as_bigmat) return(kinlist)
  
  # ----------------------------------
  # big matrix assembly by index mapping
  # ----------------------------------
  
  bigmat <- matrix(0, nsamp, nsamp)
  rownames(bigmat) <- rownames(dart_data$gt)
  colnames(bigmat) <- rownames(dart_data$gt)
  
  for (grp in names(kinlist)) {
    im <- kinlist[[grp]]
    inds <- rownames(im)
    
    # find positions in full sample order
    idx <- match(inds, rownames(bigmat))
    
    # copy sub-matrix
    bigmat[idx, idx] <- im
  }
  
  return(bigmat)
}

