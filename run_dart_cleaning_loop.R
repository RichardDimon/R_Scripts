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
    remove_pops_less_than_n5 = "FALSE", # keep your original string style handling
    samples_per_pop_remove = 5,
    downsample = "FALSE",
    samples_per_pop_downsample = 5,
    output_dir = paste0(species, "/outputs_", site_col_name, "_", species_col_name, "/"),
    max_rounds = 20 # safety cap to prevent runaway loops
) {
  
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
  
  # helper: canonicalize site names when multiple species share site name (same logic as your loop)
  disambiguate_site_names <- function(dms_obj) {
    for (g in 1:length(unique(dms_obj$meta$analyses[,site_col_name]))){
      siteval <- unique(dms_obj$meta$analyses[,site_col_name])[g]
      numz <- which(dms_obj$meta$analyses[,site_col_name] == siteval)
      if (length(unique(dms_obj$meta$analyses[,species_col_name][numz]))>1){
        for(t in 1:length(unique(dms_obj$meta$analyses[,species_col_name][numz]))){
          numz2 <- which(dms_obj$meta$analyses[,species_col_name][numz]==unique(dms_obj$meta$analyses[,species_col_name][numz])[t])
          dms_obj$meta$analyses[,site_col_name][numz][numz2] <-  paste0(dms_obj$meta$analyses[,site_col_name][numz][numz2],"_",t)
        }
      }
    }
    return(dms_obj)
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
    # fix site names where needed
    dms <- disambiguate_site_names(dms)
    
    # kinship -> identify clones (as per your code)
    kin <- individual_kinship_by_pop_sp(dart_data = dms, basedir = RandRbase, species = species, dataset = dataset, pop = dms$meta$analyses[,site_col_name], sp = dms$meta$analyses[,species_col_name], maf=maf_val, mis=locus_miss, as_bigmat=TRUE)
    kin[is.na(kin)] <- 0
    kin3 <- ifelse(kin < clonal_threshold, 0, 1)
    
    # adjacency network & clustering
    network <- igraph::graph_from_adjacency_matrix(kin3, mode="undirected", diag=F, weighted=T)
    ceb <- igraph::cluster_fast_greedy(network)
    clones <- as.data.frame(cbind(genet = ceb$membership, sample = ceb$names), stringsAsFactors = FALSE)
    clones_out <- merge(clones, dms$meta$analyses[, c("sample", "lat", "long", site_col_name, species_col_name)], by="sample")
    clones_out <- clones_out[order(as.numeric(clones_out$genet)), ]
    clones_out$genet <- as.numeric(clones_out$genet)
    
    # write clones table for this round
    openxlsx::write.xlsx(clones_out,
                         file = file.path(tables_dir, paste0("PLINK_clones_round", round, ".xlsx")),
                         asTable = FALSE, overwrite = TRUE)
    
    # identify genets with multiple ramets
    clonenum <- which(duplicated(clones_out$genet) == TRUE)
    clonenum2 <- unique(clones_out$genet[clonenum])
    
    if (length(clonenum2) == 0) {
      message("No genets with multiple ramets detected (round ", round, ").")
      log_df <- rbind(log_df, data.frame(round=round, step="clones", removed=0))
    } else {
      keepsamps <- c()
      for (w in seq_along(clonenum2)) {
        samps <- clones_out$sample[which(clones_out$genet == clonenum2[w])]
        num_missing <- data.frame(NSWNum = samps, sumNA = NA, stringsAsFactors = FALSE)
        for (s in seq_along(samps)) {
          num <- which(dms$sample_names == samps[s])
          num_missing$sumNA[s] <- sum(is.na(dms$gt[num, ]) == TRUE)
        }
        num_missing <- num_missing[order(num_missing$sumNA), ]
        keepsamps[w] <- num_missing$NSWNum[1]
      }
      
      allgenets <- 1:max(clones_out$genet)
      noclonegenets <- which(!allgenets %in% clonenum2)
      noclonesamps <- clones_out$sample[which(clones_out$genet %in% noclonegenets)]
      
      SampstoKeep <- c(keepsamps, noclonesamps) # final keep list
      dms2 <- remove.by.list(dm, SampstoKeep)
      
      removed_here <- length(dms$sample_names) - length(dms2$sample_names)
      if (removed_here < 0) removed_here <- 0 # safety
      removed_this_round_total <- removed_this_round_total + removed_here
      log_df <- rbind(log_df, data.frame(round=round, step="clones", removed=removed_here))
      message("Removed ", removed_here, " samples as clones/ramet selection (round ", round, ").")
      
      dms <- dms2
      # update treatment, site summary etc.
      treatment <- dms$treatment
      dms <- disambiguate_site_names(dms)
      unfiltered_site_summary <- dms$meta$analyses %>% as.data.frame() %>%
        group_by(!!rlang::sym(species_col_name), !!rlang::sym(site_col_name)) %>%
        dplyr::summarize(n_unfiltered = sum(n()),
                         lat = mean(as.numeric(lat), na.rm=TRUE),
                         long = mean(as.numeric(long),na.rm=TRUE),
                         .groups = 'drop') %>%
        filter(n_unfiltered > 0) %>%
        as.data.frame()
      unfiltered_site_summary <- unfiltered_site_summary[rev(order(unfiltered_site_summary$lat)),]
    }
    
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
    
    # Final clone re-check (as in your repeated script, you ran clones again)
    dms <- disambiguate_site_names(dms)
    kin_final <- individual_kinship_by_pop_sp(dms, RandRbase, species, dataset, pop = dms$meta$analyses[,site_col_name], sp = dms$meta$analyses[,species_col_name], maf=maf_val, mis=locus_miss, as_bigmat=TRUE)
    kin_final[is.na(kin_final)] <- 0
    kin3_final <- ifelse(kin_final < clonal_threshold, 0, 1)
    network_final <- igraph::graph_from_adjacency_matrix(kin3_final, mode="undirected", diag=F, weighted=T)
    ceb_final <- igraph::cluster_fast_greedy(network_final)
    clones_final <- as.data.frame(cbind(genet = ceb_final$membership, sample = ceb_final$names), stringsAsFactors = FALSE)
    clones_out_final <- merge(clones_final, dms$meta$analyses[, c("sample","lat","long",site_col_name,species_col_name)], by="sample")
    clones_out_final <- clones_out_final[order(as.numeric(clones_out_final$genet)),]
    clones_out_final$genet <- as.numeric(clones_out_final$genet)
    openxlsx::write.xlsx(clones_out_final, file = file.path(tables_dir, paste0("PLINK_clones_round", round, "_finalcheck.xlsx")), asTable = FALSE, overwrite = TRUE)
    
    clonenum_final <- which(duplicated(clones_out_final$genet) == TRUE)
    clonenum2_final <- unique(clones_out_final$genet[clonenum_final])
    if (length(clonenum2_final) > 0) {
      keepsamps_final <- c()
      for (w in seq_along(clonenum2_final)) {
        samps <- clones_out_final$sample[which(clones_out_final$genet == clonenum2_final[w])]
        num_missing <- data.frame(NSWNum = samps, sumNA = NA, stringsAsFactors = FALSE)
        for (s in seq_along(samps)) {
          num <- which(dms$sample_names == samps[s])
          num_missing$sumNA[s] <- sum(is.na(dms$gt[num, ]) == TRUE)
        }
        num_missing <- num_missing[order(num_missing$sumNA), ]
        keepsamps_final[w] <- num_missing$NSWNum[1]
      }
      allgenets <- 1:max(clones_out_final$genet)
      noclonegenets <- which(!allgenets %in% clonenum2_final)
      noclonesamps <- clones_out_final$sample[which(clones_out_final$genet %in% noclonegenets)]
      SampstoKeep_final <- c(keepsamps_final, noclonesamps)
      dms2_final <- remove.by.list(dm, SampstoKeep_final)
      
      removed_final <- length(dms$sample_names) - length(dms2_final$sample_names)
      if (removed_final < 0) removed_final <- 0
      removed_this_round_total <- removed_this_round_total + removed_final
      log_df <- rbind(log_df, data.frame(round=round, step="clones_finalcheck", removed=removed_final))
      dms <- dms2_final
      message("Final clone-check removed ", removed_final, " samples (round ", round, ").")
    } else {
      log_df <- rbind(log_df, data.frame(round=round, step="clones_finalcheck", removed=0))
      message("No clones in final check (round ", round, ").")
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
  message("Restoring original site names at final output…")
  
  # merge original site names back into final dms object
  restored <- merge(
    dms$meta$analyses,
    orig_site_lookup,
    by = "sample",
    all.x = TRUE,
    sort = FALSE
  )
  
  # replace modified disambiguated names with original
  dms$meta$analyses[, site_col_name] <- restored$site_original
  
  
  
  # final outputs
  final_list <- list(dms = dms, d1 = d1, rounds = round, log = log_df, treatment = treatment)
  return(final_list)
}



