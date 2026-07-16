# Calculate kinship separately within groups supplied through `pop`.
# Matrix assembly uses the original sample indices, so groups do not need to
# occupy consecutive rows in dart_data$gt.
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
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop("The SNPRelate package is required.", call. = FALSE)
  }

  sample_ids <- rownames(dart_data$gt)
  n_samples <- nrow(dart_data$gt)

  if (is.null(sample_ids) || any(is.na(sample_ids)) || any(sample_ids == "")) {
    stop("dart_data$gt must have non-empty sample IDs as row names.", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    stop("Sample IDs in dart_data$gt must be unique.", call. = FALSE)
  }
  if (length(pop) != n_samples) {
    stop("pop must contain one value for every genotype row.", call. = FALSE)
  }

  group_values <- as.character(pop)
  group_values[
    is.na(group_values) | trimws(group_values) == ""
  ] <- "Unassigned_population"
  group_levels <- unique(group_values)

  gds_file <- dart2gds(
    dart_data,
    basedir,
    species,
    dataset
  )

  if (!file.exists(gds_file)) {
    stop("The GDS file was not created: ", gds_file, call. = FALSE)
  }

  gds <- SNPRelate::snpgdsOpen(gds_file)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  kinlist <- vector("list", length(group_levels))
  names(kinlist) <- group_levels

  if (isTRUE(as_bigmat)) {
    output_matrix <- matrix(
      0,
      nrow = n_samples,
      ncol = n_samples,
      dimnames = list(sample_ids, sample_ids)
    )
  }

  for (group_value in group_levels) {
    group_indices <- which(group_values == group_value)
    group_samples <- sample_ids[group_indices]

    if (length(group_samples) == 1L) {
      group_matrix <- matrix(
        0.5,
        nrow = 1L,
        ncol = 1L,
        dimnames = list(group_samples, group_samples)
      )
    } else {
      group_result <- SNPRelate::snpgdsIBDMoM(
        gds,
        sample.id = group_samples,
        maf = maf,
        missing.rate = mis,
        num.thread = num.thread,
        kinship = TRUE
      )

      group_matrix <- as.matrix(group_result$kinship)
      rownames(group_matrix) <- group_samples
      colnames(group_matrix) <- group_samples
      group_matrix[is.na(group_matrix)] <- 0
      group_matrix <- pmax(group_matrix, t(group_matrix))
    }

    kinlist[[group_value]] <- group_matrix

    if (isTRUE(as_bigmat)) {
      output_matrix[group_indices, group_indices] <- group_matrix
    }
  }

  if (isTRUE(as_bigmat)) {
    return(output_matrix)
  }

  kinlist
}
# Calculate kinship using species-wide allele-frequency references, then mask
# comparisons to retain only pairs from the same species and site.
#
# This preserves a consistent SNP/MAF reference across sites within each
# species while allowing clone decisions to be restricted to within-site pairs.
individual_kinship_by_pop_sp <- function(
    dart_data,
    basedir,
    species,
    dataset,
    pop,
    sp,
    maf = 0.1,
    mis = 0.2,
    as_bigmat = TRUE,
    num.thread = 1,
    return_reference = FALSE,
    mask_value = 0
) {
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop("The SNPRelate package is required.", call. = FALSE)
  }

  sample_ids <- rownames(dart_data$gt)
  n_samples <- nrow(dart_data$gt)

  if (is.null(sample_ids) || any(is.na(sample_ids)) || any(sample_ids == "")) {
    stop("dart_data$gt must have non-empty sample IDs as row names.", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    stop("Sample IDs in dart_data$gt must be unique.", call. = FALSE)
  }
  if (length(pop) != n_samples) {
    stop("pop must contain one site value for every genotype row.", call. = FALSE)
  }
  if (length(sp) != n_samples) {
    stop("sp must contain one species value for every genotype row.", call. = FALSE)
  }
  if (!exists("individual_kinship_by_pop", mode = "function")) {
    stop(
      "individual_kinship_by_pop() must be loaded before individual_kinship_by_pop_sp().",
      call. = FALSE
    )
  }

  site_values <- as.character(pop)
  species_values <- as.character(sp)

  site_values[
    is.na(site_values) | trimws(site_values) == ""
  ] <- "Unassigned_site"
  species_values[
    is.na(species_values) | trimws(species_values) == ""
  ] <- "Unassigned_species"

  # Estimate kinship separately for each species, using all sites from that
  # species in the allele-frequency, MAF and missingness calculations.
  species_reference <- individual_kinship_by_pop(
    dart_data = dart_data,
    basedir = basedir,
    species = species,
    dataset = dataset,
    pop = species_values,
    maf = maf,
    mis = mis,
    as_bigmat = TRUE,
    num.thread = num.thread
  )

  species_reference <- as.matrix(species_reference)

  if (is.null(rownames(species_reference))) {
    rownames(species_reference) <- sample_ids
  }
  if (is.null(colnames(species_reference))) {
    colnames(species_reference) <- sample_ids
  }

  if (!all(sample_ids %in% rownames(species_reference)) ||
      !all(sample_ids %in% colnames(species_reference))) {
    stop(
      "The species-level kinship matrix does not contain all genotype samples.",
      call. = FALSE
    )
  }

  species_reference <- species_reference[
    sample_ids,
    sample_ids,
    drop = FALSE
  ]
  species_reference[is.na(species_reference)] <- 0
  species_reference <- pmax(species_reference, t(species_reference))
  dimnames(species_reference) <- list(sample_ids, sample_ids)

  # Only within-species, within-site pairs are eligible for clone decisions.
  same_species <- outer(species_values, species_values, FUN = "==")
  same_site <- outer(site_values, site_values, FUN = "==")
  comparison_mask <- same_species & same_site
  dimnames(comparison_mask) <- list(sample_ids, sample_ids)

  species_site_masked <- species_reference
  species_site_masked[!comparison_mask] <- mask_value
  dimnames(species_site_masked) <- list(sample_ids, sample_ids)

  # Confirm that masking has not altered any eligible within-site estimate.
  eligible_difference <- species_site_masked[comparison_mask] -
    species_reference[comparison_mask]

  if (length(eligible_difference) > 0L &&
      any(abs(eligible_difference) > sqrt(.Machine$double.eps), na.rm = TRUE)) {
    stop(
      "Masking unexpectedly changed one or more within-species, within-site kinship values.",
      call. = FALSE
    )
  }

  if (isTRUE(return_reference)) {
    return(list(
      species_site = species_site_masked,
      species_reference = species_reference,
      comparison_mask = comparison_mask,
      species = species_values,
      site = site_values
    ))
  }

  if (isTRUE(as_bigmat)) {
    return(species_site_masked)
  }

  # Preserve the historical list-style output when as_bigmat = FALSE, but
  # values are now drawn from the species-wide reference matrix.
  species_site_groups <- paste(species_values, site_values, sep = "__")
  group_levels <- unique(species_site_groups)

  kinlist <- lapply(group_levels, function(group_value) {
    idx <- which(species_site_groups == group_value)
    species_site_masked[idx, idx, drop = FALSE]
  })
  names(kinlist) <- group_levels

  kinlist
}
