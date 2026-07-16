run_dart_cleaning_loop <- function(
    d1,                     # initial DArT object (your d1)
    dms,                    # initial dms (dms2 from your setup)
    RandRbase,
    species,
    dataset,
    species_col_name,
    site_col_name,
    sample_miss,
    locus_miss,
    maf_val,
    clonal_threshold,
    clone_scope = c("site", "species"), # "site" = within species × site; "species" = across sites within species
    remove_pops_less_than_n5 = "FALSE", # keep your original string style handling
    samples_per_pop_remove = 5,
    downsample = "FALSE",
    samples_per_pop_downsample = 5,
    output_dir = paste0(species, "/outputs_", site_col_name, "_", species_col_name, "/"),
    max_rounds = 20 # safety cap to prevent runaway loops
) {
  
  clone_scope <- tolower(trimws(as.character(clone_scope[1])))
  if (!(clone_scope %in% c("site", "species"))) {
    stop("clone_scope must be either 'site' or 'species'.", call. = FALSE)
  }
  message("Clone-removal scope: ", clone_scope)
  
  orig_site_names <- dms$meta$analyses[, site_col_name]
  orig_site_lookup <- data.frame(
    sample = dms$sample_names,
    site_original = orig_site_names,
    stringsAsFactors = FALSE
  )
  
  
  
  # helper: ensure output subdirs exist
  plots_dir  <- file.path(output_dir, "plots")
  tables_dir <- file.path(output_dir, "tables")
  rfiles_dir <- file.path(output_dir, "r_files")
  for (p in c(plots_dir, tables_dir, rfiles_dir)){
    if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  }
  
  round <- 0
  prev_n <- length(dms$sample_names)
  log_df <- data.frame(round = integer(),
                       step = character(),
                       removed = integer(),
                       stringsAsFactors = FALSE)
  
  # small helper to run reanalysis (d1 -> d2 -> d3 -> m1 -> dm -> dmsnew)
  reanalyse <- function(d1_current, dms_current, round_tag) {
    # remove from d1 any samples present in dms_current (this matches your remove.by.list usage)
    d1old <- d1_current
    d1_new <- remove.by.list(d1_current, dms_current$sample_names)
    qc1 <- report.dart.qc.stats.RD(d1_new, RandRbase, species, dataset, threshold_missing_loci = sample_miss)
    d2_new <- remove.poor.quality.snps(d1_new, min_repro=0.96, max_missing=locus_miss)
    qc2 <- report.dart.qc.stats.RD(d2_new, RandRbase, species, dataset, threshold_missing_loci = sample_miss)
    d3_new <- sample.one.snp.per.locus.random(d2_new, seed=214241)
    d3_new$treatment <- paste0(d3_new$treatment, "_rep", round_tag)
    qc3 <- report.dart.qc.stats.RD(d3_new, RandRbase, species, dataset, threshold_missing_loci = sample_miss)
    
    # meta
    m1 <- read.meta.data.full.analyses.df(d3_new, RandRbase, species, dataset)
    m1$analyses[,site_col_name] <- gsub(" ", "_", m1$analyses[,site_col_name])
    orderbylat <- rev(order(m1$lat))
    m1$lat <- m1$lat[orderbylat]
    m1$sample_names <- m1$sample_names[orderbylat]
    m1$site <- m1$site[orderbylat]
    m1$long <- m1$long[orderbylat]
    m1$analyses <- m1$analyses[orderbylat,]
    # write any missing meta matches (as in your original code)
    write.table(data.frame(sample=d3_new$sample_names[!(d3_new$sample_names %in% m1$sample_names)]),
                file = file.path(tables_dir, "samples_in_dart_not_in_meta.tsv"),
                sep = "\t", row.names = FALSE)
    
    dm_new <- dart.meta.data.merge(d3_new, m1)
    
    # remove samples that are NA for site variable
    samples_with_site_variable <- dm_new$sample_names[!is.na(dm_new$meta$analyses[,site_col_name])]
    dmsnew <- remove.by.list(dm_new, samples_with_site_variable)
    
    return(list(d1=d1_new, d2=d2_new, d3=d3_new, dm=dm_new, dms=dmsnew,
                qc1=qc1, qc2=qc2, qc3=qc3))
  }
  
  # helper: canonicalize site names when the same site label occurs in multiple species
  # This is retained for downstream site-based summaries and filtering.
  disambiguate_site_names <- function(dms_obj) {
    meta <- as.data.frame(dms_obj$meta$analyses, stringsAsFactors = FALSE)
    site_values <- as.character(meta[[site_col_name]])
    species_values <- as.character(meta[[species_col_name]])

    for (siteval in unique(site_values)) {
      if (is.na(siteval)) next
      site_idx <- which(site_values == siteval)
      site_species <- unique(species_values[site_idx])
      site_species <- site_species[!is.na(site_species)]

      if (length(site_species) > 1) {
        for (i in seq_along(site_species)) {
          target_idx <- site_idx[species_values[site_idx] == site_species[i]]
          site_values[target_idx] <- paste0(siteval, "_", i)
        }
      }
    }

    dms_obj$meta$analyses[, site_col_name] <- site_values
    dms_obj
  }

  # helper: calculate a kinship matrix using the requested clone-removal scope
  calculate_clone_kinship <- function(dms_obj) {
    meta <- as.data.frame(dms_obj$meta$analyses, stringsAsFactors = FALSE)

    if (!"sample" %in% names(meta)) {
      meta$sample <- dms_obj$sample_names
    }

    meta_match <- match(dms_obj$sample_names, meta$sample)
    if (anyNA(meta_match)) {
      stop(
        "Some dms samples could not be matched to dms$meta$analyses$sample during clone detection.",
        call. = FALSE
      )
    }
    meta <- meta[meta_match, , drop = FALSE]

    species_values <- as.character(meta[[species_col_name]])
    site_values <- as.character(meta[[site_col_name]])

    species_values[is.na(species_values) | species_values == ""] <- "Unassigned_species"
    site_values[is.na(site_values) | site_values == ""] <- "Unassigned_site"

    if (clone_scope == "site") {
      # Site scope: compare samples only within each species × site combination.
      if (!exists("individual_kinship_by_pop_sp", mode = "function")) {
        stop(
          "clone_scope = 'site' requires individual_kinship_by_pop_sp().",
          call. = FALSE
        )
      }

      kin <- individual_kinship_by_pop_sp(
        dart_data = dms_obj,
        basedir = RandRbase,
        species = species,
        dataset = dataset,
        pop = site_values,
        sp = species_values,
        maf = maf_val,
        mis = locus_miss,
        as_bigmat = TRUE
      )
    } else {
      # Species scope: compare samples across all sites, but only within species.
      if (!exists("individual_kinship_by_pop", mode = "function")) {
        stop(
          "clone_scope = 'species' requires individual_kinship_by_pop().",
          call. = FALSE
        )
      }

      kin <- individual_kinship_by_pop(
        dms_obj,
        RandRbase,
        species,
        dataset,
        species_values,
        maf = maf_val,
        mis = locus_miss,
        as_bigmat = TRUE
      )
    }

    kin <- as.matrix(kin)
    sample_names <- dms_obj$sample_names

    if (is.null(rownames(kin))) rownames(kin) <- sample_names
    if (is.null(colnames(kin))) colnames(kin) <- sample_names

    if (!all(sample_names %in% rownames(kin)) ||
        !all(sample_names %in% colnames(kin))) {
      stop(
        "The kinship matrix sample names do not match dms$sample_names.",
        call. = FALSE
      )
    }

    kin <- kin[sample_names, sample_names, drop = FALSE]
    kin[is.na(kin)] <- 0

    # Protect against small numerical asymmetries in the returned matrix.
    kin <- pmax(kin, t(kin))
    dimnames(kin) <- list(sample_names, sample_names)
    diag(kin) <- 0
    kin
  }

  # helper: identify connected clone groups and retain the sample with the fewest missing loci
  filter_clones_once <- function(dms_obj, removal_source_obj, round_number, stage = "main") {
    if (length(dms_obj$sample_names) <= 1) {
      return(list(
        dms = dms_obj,
        removed = 0L,
        removed_samples = character(),
        assignments = data.frame(),
        kinship = NULL
      ))
    }

    kin <- calculate_clone_kinship(dms_obj)

    clone_adj <- kin >= clonal_threshold
    clone_adj[is.na(clone_adj)] <- FALSE
    diag(clone_adj) <- FALSE

    network <- igraph::graph_from_adjacency_matrix(
      clone_adj * 1L,
      mode = "undirected",
      diag = FALSE,
      weighted = NULL
    )

    clone_membership <- igraph::components(network)$membership
    clones <- data.frame(
      sample = names(clone_membership),
      genet = as.integer(unname(clone_membership)),
      stringsAsFactors = FALSE
    )

    meta <- as.data.frame(dms_obj$meta$analyses, stringsAsFactors = FALSE)
    if (!"sample" %in% names(meta)) meta$sample <- dms_obj$sample_names
    meta <- meta[match(clones$sample, meta$sample), , drop = FALSE]

    missing_counts <- rowSums(is.na(dms_obj$gt))
    names(missing_counts) <- dms_obj$sample_names

    clones_out <- data.frame(
      genet = clones$genet,
      sample = clones$sample,
      n_missing_loci = as.integer(missing_counts[clones$sample]),
      stringsAsFactors = FALSE
    )

    metadata_columns <- intersect(
      c("lat", "long", site_col_name, species_col_name),
      names(meta)
    )
    if (length(metadata_columns) > 0) {
      clones_out <- cbind(
        clones_out,
        meta[, metadata_columns, drop = FALSE]
      )
    }

    group_sizes <- table(clones_out$genet)
    clones_out$group_size <- as.integer(group_sizes[as.character(clones_out$genet)])
    clones_out$clone_detected <- clones_out$group_size > 1
    clones_out$retained <- TRUE

    clone_genets <- as.integer(names(group_sizes[group_sizes > 1]))

    for (genet_id in clone_genets) {
      group_idx <- which(clones_out$genet == genet_id)
      group_order <- order(
        clones_out$n_missing_loci[group_idx],
        clones_out$sample[group_idx]
      )
      retained_idx <- group_idx[group_order[1]]
      clones_out$retained[group_idx] <- FALSE
      clones_out$retained[retained_idx] <- TRUE
    }

    clones_out$clone_scope <- clone_scope
    clones_out$clone_group <- if (clone_scope == "site") {
      paste(
        clones_out[[species_col_name]],
        clones_out[[site_col_name]],
        sep = "__"
      )
    } else {
      as.character(clones_out[[species_col_name]])
    }
    clones_out$removal_reason <- ifelse(
      clones_out$clone_detected & !clones_out$retained,
      paste0("clone_removed_", clone_scope),
      ifelse(clones_out$clone_detected, "clone_representative_retained", "not_a_clone")
    )

    clones_out <- clones_out[order(clones_out$genet, !clones_out$retained,
                                   clones_out$n_missing_loci, clones_out$sample), ]

    stage_suffix <- if (identical(stage, "main")) "" else paste0("_", stage)
    output_file <- file.path(
      tables_dir,
      paste0(
        "PLINK_clones_", clone_scope,
        "_round", round_number,
        stage_suffix,
        ".xlsx"
      )
    )

    openxlsx::write.xlsx(
      clones_out,
      file = output_file,
      asTable = FALSE,
      overwrite = TRUE
    )

    samples_to_keep <- clones_out$sample[clones_out$retained]
    removed_samples <- clones_out$sample[!clones_out$retained]

    filtered_obj <- remove.by.list(removal_source_obj, samples_to_keep)
    removed_n <- length(dms_obj$sample_names) - length(filtered_obj$sample_names)

    if (removed_n < 0) {
      stop(
        "Clone filtering increased the number of samples, indicating unexpected remove.by.list behaviour.",
        call. = FALSE
      )
    }

    if (removed_n != length(removed_samples)) {
      warning(
        "The number removed by remove.by.list (", removed_n,
        ") differs from the clone-decision table (", length(removed_samples), ")."
      )
    }

    list(
      dms = filtered_obj,
      removed = as.integer(removed_n),
      removed_samples = removed_samples,
      assignments = clones_out,
      kinship = kin
    )
  }

  # main loop
  while (TRUE) {
    round <- round + 1
    if (round > max_rounds) {
      message("Reached max_rounds cap (", max_rounds, "). Stopping loop.")
      break
    }
    message("=== ROUND ", round, " ===")
    removed_this_round_total <- 0
    
    # -------------------------
    # 1) Remove high missingness samples
    # -------------------------
    samples_high_missing <- dms$sample_names[which(rowMeans(is.na(dms$gt)) > sample_miss)]
    if (length(samples_high_missing) > 0) {
      # write table with round number
      write.table(data.frame(sample=samples_high_missing,
                             sample_miss = rowMeans(is.na(dms$gt))[which(rowMeans(is.na(dms$gt))>sample_miss)]),
                  file = file.path(tables_dir, paste0("high_missing_samples_removed_round", round, ".tsv")),
                  sep = "\t", row.names = FALSE)
      dms <- remove.by.missingness(dms, sample_miss)
      removed_this_round_total <- removed_this_round_total + length(samples_high_missing)
      log_df <- rbind(log_df, data.frame(round=round, step="missingness", removed=length(samples_high_missing)))
      message("Removed ", length(samples_high_missing), " samples for high missingness (round ", round, ").")
    } else {
      message("No high-missingness samples found (round ", round, ").")
      log_df <- rbind(log_df, data.frame(round=round, step="missingness", removed=0))
    }
    
    # Reanalyse after missingness (use round number as tag)
    re_res <- reanalyse(d1, dms, paste0(round, "a"))
    d1 <- re_res$d1
    d2 <- re_res$d2
    d3 <- re_res$d3
    dm <- re_res$dm
    dms <- re_res$dms
    
    # compute unfiltered_site_summary used later
    unfiltered_site_summary <- dms$meta$analyses %>% as.data.frame() %>%
      group_by(!!rlang::sym(species_col_name), !!rlang::sym(site_col_name)) %>%
      dplyr::summarize(n_unfiltered = sum(n()),
                       lat = mean(as.numeric(lat), na.rm=TRUE),
                       long = mean(as.numeric(long), na.rm=TRUE),
                       .groups = 'drop') %>%
      filter(n_unfiltered > 0) %>%
      as.data.frame()
    unfiltered_site_summary <- unfiltered_site_summary[rev(order(unfiltered_site_summary$lat)),]
    
    # -------------------------
    # 2) Clone detection & removal
    # -------------------------
    dms <- disambiguate_site_names(dms)

    clone_result <- filter_clones_once(
      dms_obj = dms,
      removal_source_obj = dm,
      round_number = round,
      stage = "main"
    )

    dms <- clone_result$dms
    removed_here <- clone_result$removed
    removed_this_round_total <- removed_this_round_total + removed_here

    log_df <- rbind(
      log_df,
      data.frame(
        round = round,
        step = paste0("clones_", clone_scope),
        removed = removed_here,
        stringsAsFactors = FALSE
      )
    )

    if (removed_here > 0) {
      message(
        "Removed ", removed_here,
        " clone samples using the '", clone_scope,
        "' scope (round ", round, ")."
      )
    } else {
      message(
        "No clone samples detected using the '", clone_scope,
        "' scope (round ", round, ")."
      )
    }

    treatment <- dms$treatment
    dms <- disambiguate_site_names(dms)

    unfiltered_site_summary <- dms$meta$analyses %>%
      as.data.frame() %>%
      group_by(
        !!rlang::sym(species_col_name),
        !!rlang::sym(site_col_name)
      ) %>%
      dplyr::summarize(
        n_unfiltered = dplyr::n(),
        lat = mean(as.numeric(lat), na.rm = TRUE),
        long = mean(as.numeric(long), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_unfiltered > 0) %>%
      as.data.frame()

    unfiltered_site_summary <- unfiltered_site_summary[
      rev(order(unfiltered_site_summary$lat)),
    ]

    # Reanalyse after clone removal
    re_res2 <- reanalyse(d1, dms, paste0(round, "b"))
    d1 <- re_res2$d1
    d2 <- re_res2$d2
    d3 <- re_res2$d3
    dm <- re_res2$dm
    dms <- re_res2$dms
    
    # -------------------------
    # 3) Pop filtering and downsampling
    # -------------------------
    if (toupper(remove_pops_less_than_n5) == "TRUE") {
      not_n5_sites <- as.vector(names(which(table(dms$meta$analyses[, site_col_name]) < samples_per_pop_remove)))
      not_n5_samples <- dms$sample_names[which(!(dms$meta$analyses[, site_col_name] %in% not_n5_sites))]
      dms <- remove.by.list(dms, not_n5_samples)
      log_df <- rbind(log_df, data.frame(round=round, step="remove_small_pops", removed=length(not_n5_sites)))
      message("Removed sites with < ", samples_per_pop_remove, " (", length(not_n5_sites), " sites) (round ", round, ").")
    } else {
      # create dms_no_n1_sites as your code used it (no removal from dms in this branch)
      not_n1_sites <- as.vector(unfiltered_site_summary[unfiltered_site_summary$n_unfiltered <= 1, 2])
      not_n1_samples <- dms$sample_names[which(!(dms$meta$analyses[, site_col_name] %in% not_n1_sites))]
      dms_no_n1_sites <- remove.by.list(dms, not_n1_samples)
      # we don't assign back to dms unless user requested removal above
      log_df <- rbind(log_df, data.frame(round=round, step="remove_small_pops", removed=0))
    }
    
    
    
    
    downsample_max_diversity <- function(geno_mat, sample_ids, n_keep) {
      # geno_mat: numeric genotype matrix individuals × SNPs
      # sample_ids: vector of sample names (same order as rows in geno_mat)
      # n_keep: number of individuals to retain
      
      if(length(sample_ids) <= n_keep) {
        return(sample_ids)  # nothing to do
      }
      
      # Compute Euclidean genetic distances
      d <- as.matrix(dist(geno_mat))
      diag(d) <- NA
      
      # 1. Pick the individual with the maximum average distance to others
      avg_dist <- rowMeans(d, na.rm = TRUE)
      first <- sample_ids[which.max(avg_dist)]
      selected <- c(first)
      
      # 2. Iteratively add the individual furthest from the selected set
      remaining <- setdiff(sample_ids, selected)
      
      while(length(selected) < n_keep) {
        # For each remaining individual, compute its distance to the nearest selected one
        min_dists <- sapply(remaining, function(ind) {
          inds <- which(sample_ids %in% selected)
          rem  <- which(sample_ids == ind)
          min(d[rem, inds], na.rm = TRUE)
        })
        
        # Pick the individual with the largest minimum distance
        next_ind <- remaining[which.max(min_dists)]
        selected <- c(selected, next_ind)
        remaining <- setdiff(remaining, next_ind)
      }
      return(selected)
    }
    
    
    if (toupper(downsample) == "TRUE") {
      geno <- dms$gt  # or however your genotype matrix is stored
      samples_to_keep <- c()
      sites <- unique(dms$meta$analyses[, site_col_name])
      
      for (s in sites) {
        site_samples <- dms$sample_names[dms$meta$analyses[, site_col_name] == s]
        if (length(site_samples) > samples_per_pop_downsample) {
          # extract genotype rows
          geno_sub <- geno[site_samples, , drop = FALSE]
          keep <- downsample_max_diversity(geno_mat = geno_sub, sample_ids = site_samples, n_keep = samples_per_pop_downsample)
          samples_to_keep <- c(samples_to_keep, keep)
        } else {
          samples_to_keep <- c(samples_to_keep, site_samples)
        }
      }
      # remove all samples not in the genetically diverse kept set
      dms <- remove.by.list(dms, samples_to_keep)
      message("Downsampled populations to ", samples_per_pop_downsample, " (round ", round, ").")
    } else {
      log_df <- rbind(log_df, data.frame(round=round, step="downsample", removed=0))
    }
    
    # Reanalyse after pop filtering / downsample
    re_res3 <- reanalyse(d1, dms, paste0(round, "c"))
    d1 <- re_res3$d1
    d2 <- re_res3$d2
    d3 <- re_res3$d3
    dm <- re_res3$dm
    dms <- re_res3$dms
    
    # Final clone re-check using the same selected scope
    dms <- disambiguate_site_names(dms)

    clone_final_result <- filter_clones_once(
      dms_obj = dms,
      removal_source_obj = dm,
      round_number = round,
      stage = "finalcheck"
    )

    dms <- clone_final_result$dms
    removed_final <- clone_final_result$removed
    removed_this_round_total <- removed_this_round_total + removed_final

    log_df <- rbind(
      log_df,
      data.frame(
        round = round,
        step = paste0("clones_", clone_scope, "_finalcheck"),
        removed = removed_final,
        stringsAsFactors = FALSE
      )
    )

    if (removed_final > 0) {
      message(
        "Final clone check removed ", removed_final,
        " samples using the '", clone_scope,
        "' scope (round ", round, ")."
      )
    } else {
      message(
        "No clones detected in the final '", clone_scope,
        "' check (round ", round, ")."
      )
    }

    # Write a small round summary file
    round_summary <- data.frame(
      round = round,
      samples_after = length(dms$sample_names),
      removed_this_round = removed_this_round_total
    )
    write.table(round_summary, file = file.path(tables_dir, paste0("round_summary_", round, ".tsv")), sep = "\t", row.names = FALSE)
    
    # stopping condition: stop when sample count doesn't change
    new_n <- length(dms$sample_names)
    message("Round ", round, " end: samples = ", new_n, " (previous ", prev_n, ").")
    if (new_n == prev_n) {
      message("No change in sample count; stopping.")
      break
    } else {
      prev_n <- new_n
      # continue looping
    }
  } # end while
  
  
  
  
  # ==================================
  # RESTORE ORIGINAL SITE NAMES
  # ==================================
  message("Restoring original site names at final output...")

  final_meta <- as.data.frame(dms$meta$analyses, stringsAsFactors = FALSE)
  if (!"sample" %in% names(final_meta)) final_meta$sample <- dms$sample_names

  site_match <- match(final_meta$sample, orig_site_lookup$sample)
  if (anyNA(site_match)) {
    warning(
      "Some final samples could not be matched to their original site names: ",
      paste(final_meta$sample[is.na(site_match)], collapse = ", ")
    )
  }

  dms$meta$analyses[, site_col_name] <-
    orig_site_lookup$site_original[site_match]

  # final outputs
  final_list <- list(
    dms = dms,
    d1 = d1,
    rounds = round,
    log = log_df,
    treatment = treatment,
    clone_scope = clone_scope
  )
  return(final_list)
}
