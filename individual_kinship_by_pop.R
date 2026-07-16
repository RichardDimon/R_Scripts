individual_kinship_by_pop <- function(
    dart_data,
    basedir,
    species,
    dataset,
    pop,
    maf = 0.1,
    mis = 0.2,
    as_bigmat = TRUE,
    num.thread = 1
) {
  require(SNPRelate)
  
  sample_ids <- rownames(dart_data$gt)
  nsamp <- nrow(dart_data$gt)
  
  if (is.null(sample_ids) || any(sample_ids == "")) {
    stop("dart_data$gt must have non-empty sample names as row names.", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    stop("Sample names in dart_data$gt must be unique.", call. = FALSE)
  }
  if (length(pop) != nsamp) {
    stop("pop must contain one value for every row of dart_data$gt.", call. = FALSE)
  }
  
  pop <- as.character(pop)
  pop[is.na(pop) | pop == ""] <- "Unassigned_population"
  pop_levels <- unique(pop)
  
  igds_file <- dart2gds(dart_data, basedir, species, dataset)
  igds <- snpgdsOpen(igds_file)
  on.exit(snpgdsClose(igds), add = TRUE)
  
  kinlist <- vector("list", length(pop_levels))
  names(kinlist) <- pop_levels
  
  if (isTRUE(as_bigmat)) {
    bigmat <- matrix(
      0,
      nrow = nsamp,
      ncol = nsamp,
      dimnames = list(sample_ids, sample_ids)
    )
  }
  
  for (grp in pop_levels) {
    idx <- which(pop == grp)
    group_samples <- sample_ids[idx]
    
    if (length(idx) == 1L) {
      # No pairwise relationship can be estimated for a singleton group.
      ikout <- matrix(
        0,
        nrow = 1L,
        ncol = 1L,
        dimnames = list(group_samples, group_samples)
      )
    } else {
      iout <- snpgdsIBDMoM(
        igds,
        sample.id = group_samples,
        maf = maf,
        missing.rate = mis,
        num.thread = num.thread,
        kinship = TRUE
      )
      
      ikout <- iout$kinship
      rownames(ikout) <- group_samples
      colnames(ikout) <- group_samples
    }
    
    kinlist[[grp]] <- ikout
    
    if (isTRUE(as_bigmat)) {
      # Insert by the original sample indices rather than assuming that
      # population groups occupy consecutive rows in dart_data$gt.
      bigmat[idx, idx] <- ikout
    }
  }
  
  if (isTRUE(as_bigmat)) {
    return(bigmat)
  }
  
  kinlist
}