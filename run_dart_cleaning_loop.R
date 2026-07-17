run_dart_cleaning_loop <- function(
    d1,                     # initial DArT object
    dms,                    # initial merged and site-filtered DArT object
    RandRbase,
    species,
    dataset,
    species_col_name,
    site_col_name,
    sample_miss,
    locus_miss,
    maf_val,
    clonal_threshold,
    clone_scope = c("site", "species"),
    remove_pops_less_than_n5 = "FALSE",
    samples_per_pop_remove = 5,
    downsample = "FALSE",
    samples_per_pop_downsample = 5,
    output_dir = paste0(
      species, "/outputs_", site_col_name, "_", species_col_name, "/"
    ),
    max_rounds = 20,
    write_intermediate_qc = FALSE,
    write_clone_tables = TRUE,
    clone_table_include_singletons = TRUE,
    write_round_summaries = TRUE,
    save_kinship_matrices = TRUE,
    save_reference_kinship_matrices = TRUE,
    write_kinship_heatmaps = TRUE,
    kinship_show_names_max_n = 250,
    kinship_png_width_px = 1600,
    kinship_png_height_px = 1200,
    kinship_png_res = 96,
    compress_kinship_matrices = FALSE
) {

  # ============================================================
  # Validate options
  # ============================================================
  clone_scope <- tolower(trimws(as.character(clone_scope[1])))
  if (!(clone_scope %in% c("site", "species"))) {
    stop("clone_scope must be either 'site' or 'species'.", call. = FALSE)
  }

  remove_pops_flag <- identical(
    toupper(trimws(as.character(remove_pops_less_than_n5[1]))),
    "TRUE"
  )
  downsample_flag <- identical(
    toupper(trimws(as.character(downsample[1]))),
    "TRUE"
  )

  message(
    "Clone-removal scope: ",
    if (clone_scope == "site") {
      paste0(
        "site (species-wide kinship reference, with clone comparisons ",
        "restricted to the same species and site)"
      )
    } else {
      "species (species-wide kinship and clone comparisons)"
    }
  )
  message("Intermediate QC reports: ", write_intermediate_qc)

  if (!is.numeric(kinship_png_width_px) ||
      length(kinship_png_width_px) != 1L ||
      !is.finite(kinship_png_width_px) ||
      kinship_png_width_px < 400) {
    stop("kinship_png_width_px must be a single finite value >= 400.", call. = FALSE)
  }

  if (!is.numeric(kinship_png_height_px) ||
      length(kinship_png_height_px) != 1L ||
      !is.finite(kinship_png_height_px) ||
      kinship_png_height_px < 400) {
    stop("kinship_png_height_px must be a single finite value >= 400.", call. = FALSE)
  }

  if (!is.numeric(kinship_png_res) ||
      length(kinship_png_res) != 1L ||
      !is.finite(kinship_png_res) ||
      kinship_png_res <= 0) {
    stop("kinship_png_res must be a single positive finite value.", call. = FALSE)
  }

  kinship_png_width_px <- as.integer(round(kinship_png_width_px))
  kinship_png_height_px <- as.integer(round(kinship_png_height_px))
  kinship_png_res <- as.integer(round(kinship_png_res))

  # ============================================================
  # Preserve original site names for the final object
  # ============================================================
  initial_meta <- as.data.frame(dms$meta$analyses, stringsAsFactors = FALSE)
  if (!"sample" %in% names(initial_meta)) {
    initial_meta$sample <- dms$sample_names
  }

  initial_match <- match(dms$sample_names, initial_meta$sample)
  if (anyNA(initial_match)) {
    stop(
      "Some starting dms samples could not be matched to metadata.",
      call. = FALSE
    )
  }
  initial_meta <- initial_meta[initial_match, , drop = FALSE]

  orig_site_lookup <- data.frame(
    sample = dms$sample_names,
    site_original = initial_meta[[site_col_name]],
    stringsAsFactors = FALSE
  )

  # ============================================================
  # Output folders
  # ============================================================
  plots_dir          <- file.path(output_dir, "plots")
  tables_dir         <- file.path(output_dir, "tables")
  rfiles_dir         <- file.path(output_dir, "r_files")
  kinship_dir        <- file.path(rfiles_dir, "kinship_matrices")
  kinship_plots_dir  <- file.path(plots_dir, "kinship_matrices")

  output_paths <- c(plots_dir, tables_dir, rfiles_dir)
  if (isTRUE(save_kinship_matrices)) {
    output_paths <- c(output_paths, kinship_dir)
  }
  if (isTRUE(write_kinship_heatmaps)) {
    output_paths <- c(output_paths, kinship_plots_dir)
  }

  for (p in unique(output_paths)) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  }

  if (isTRUE(write_kinship_heatmaps)) {
    required_heatmap_packages <- c("ComplexHeatmap", "circlize", "colorspace")
    missing_heatmap_packages <- required_heatmap_packages[
      !vapply(required_heatmap_packages, requireNamespace, logical(1), quietly = TRUE)
    ]

    if (length(missing_heatmap_packages) > 0L) {
      stop(
        "Kinship heatmaps require these packages: ",
        paste(missing_heatmap_packages, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # ============================================================
  # State
  # ============================================================
  round <- 0L
  prev_n <- length(dms$sample_names)
  treatment <- dms$treatment

  # A final clone removal occurs after the last reanalysis in a round.
  # This flag ensures that the next round reanalyses the reduced sample set.
  needs_reanalysis <- FALSE

  log_rows <- list()
  log_index <- 0L
  kinship_files <- character()
  kinship_plot_files <- character()

  # Kinship calculations may occur several times during filtering, but only
  # the initial and final matrices/plots are retained.
  initial_kinship_saved <- FALSE
  final_kinship_candidate_ready <- FALSE
  final_candidate_samples <- NULL
  final_candidate_loci <- NULL

  initial_decision_file <- file.path(
    kinship_dir,
    paste0("kinship_", clone_scope, "_initial.rds")
  )
  final_decision_file <- file.path(
    kinship_dir,
    paste0("kinship_", clone_scope, "_final.rds")
  )
  final_reference_file <- file.path(
    kinship_dir,
    "kinship_species_reference_final.rds"
  )

  # When matrices are not being retained but a final heatmap is requested,
  # use a temporary RDS file so the final matrix does not need to remain in
  # memory throughout the loop.
  final_candidate_cache_file <- if (isTRUE(save_kinship_matrices)) {
    final_decision_file
  } else {
    tempfile(
      pattern = paste0("kinship_", clone_scope, "_final_"),
      fileext = ".rds"
    )
  }

  add_log <- function(round_number, step_name, removed_number) {
    log_index <<- log_index + 1L
    log_rows[[log_index]] <<- data.frame(
      round = as.integer(round_number),
      step = as.character(step_name),
      removed = as.integer(removed_number),
      stringsAsFactors = FALSE
    )
  }

  # ============================================================
  # Rebuild loci and metadata after the sample set changes
  # ============================================================
  reanalyse <- function(d1_current, dms_current, round_tag) {

    # remove.by.list() is used as a keep-list function in this workflow.
    d1_new <- remove.by.list(d1_current, dms_current$sample_names)

    if (isTRUE(write_intermediate_qc)) {
      qc1 <- report.dart.qc.stats.RD(
        d1_new,
        RandRbase,
        species,
        dataset,
        threshold_missing_loci = sample_miss
      )
    } else {
      qc1 <- NULL
    }

    d2_new <- remove.poor.quality.snps(
      d1_new,
      min_repro = 0.96,
      max_missing = locus_miss
    )

    if (isTRUE(write_intermediate_qc)) {
      qc2 <- report.dart.qc.stats.RD(
        d2_new,
        RandRbase,
        species,
        dataset,
        threshold_missing_loci = sample_miss
      )
    } else {
      qc2 <- NULL
    }

    d3_new <- sample.one.snp.per.locus.random(d2_new, seed = 214241)
    d3_new$treatment <- paste0(d3_new$treatment, "_rep", round_tag)

    if (isTRUE(write_intermediate_qc)) {
      qc3 <- report.dart.qc.stats.RD(
        d3_new,
        RandRbase,
        species,
        dataset,
        threshold_missing_loci = sample_miss
      )
    } else {
      qc3 <- NULL
    }

    # Read and attach metadata only when a reanalysis is genuinely required.
    m1 <- read.meta.data.full.analyses.df(
      d3_new,
      RandRbase,
      species,
      dataset
    )

    m1$analyses[, site_col_name] <- gsub(
      " ",
      "_",
      m1$analyses[, site_col_name]
    )

    orderbylat <- rev(order(m1$lat))
    m1$lat <- m1$lat[orderbylat]
    m1$sample_names <- m1$sample_names[orderbylat]
    m1$site <- m1$site[orderbylat]
    m1$long <- m1$long[orderbylat]
    m1$analyses <- m1$analyses[orderbylat, , drop = FALSE]

    missing_meta_samples <- d3_new$sample_names[
      !(d3_new$sample_names %in% m1$sample_names)
    ]

    # Avoid repeatedly writing an empty file.
    if (length(missing_meta_samples) > 0L) {
      write.table(
        data.frame(
          sample = missing_meta_samples,
          stringsAsFactors = FALSE
        ),
        file = file.path(tables_dir, "samples_in_dart_not_in_meta.tsv"),
        sep = "\t",
        row.names = FALSE
      )
    }

    dm_new <- dart.meta.data.merge(d3_new, m1)

    samples_with_site_variable <- dm_new$sample_names[
      !is.na(dm_new$meta$analyses[, site_col_name])
    ]
    dms_new <- remove.by.list(dm_new, samples_with_site_variable)

    list(
      d1 = d1_new,
      d2 = d2_new,
      d3 = d3_new,
      dm = dm_new,
      dms = dms_new,
      qc1 = qc1,
      qc2 = qc2,
      qc3 = qc3
    )
  }

  # ============================================================
  # Make site labels unique where one label occurs in >1 species
  # ============================================================
  disambiguate_site_names <- function(dms_obj) {
    meta <- as.data.frame(dms_obj$meta$analyses, stringsAsFactors = FALSE)
    site_values <- as.character(meta[[site_col_name]])
    species_values <- as.character(meta[[species_col_name]])

    for (siteval in unique(site_values)) {
      if (is.na(siteval)) next

      site_idx <- which(site_values == siteval)
      site_species <- unique(species_values[site_idx])
      site_species <- site_species[!is.na(site_species)]

      if (length(site_species) > 1L) {
        for (i in seq_along(site_species)) {
          target_idx <- site_idx[
            species_values[site_idx] == site_species[i]
          ]
          site_values[target_idx] <- paste0(siteval, "_", i)
        }
      }
    }

    dms_obj$meta$analyses[, site_col_name] <- site_values
    dms_obj
  }

  # ============================================================
  # Calculate species-referenced kinship for clone detection
  # ============================================================
  calculate_clone_kinship <- function(dms_obj) {
    meta <- as.data.frame(dms_obj$meta$analyses, stringsAsFactors = FALSE)
    if (!"sample" %in% names(meta)) meta$sample <- dms_obj$sample_names

    meta_match <- match(dms_obj$sample_names, meta$sample)
    if (anyNA(meta_match)) {
      stop(
        paste0(
          "Some dms samples could not be matched to metadata during ",
          "clone detection."
        ),
        call. = FALSE
      )
    }
    meta <- meta[meta_match, , drop = FALSE]

    species_values <- as.character(meta[[species_col_name]])

    # Recover original site values by sample name. Site labels may have been
    # temporarily disambiguated elsewhere in the workflow.
    original_site_match <- match(
      dms_obj$sample_names,
      orig_site_lookup$sample
    )
    if (anyNA(original_site_match)) {
      stop(
        paste0(
          "Some samples could not be matched to their original site values ",
          "during clone detection."
        ),
        call. = FALSE
      )
    }

    site_values <- as.character(
      orig_site_lookup$site_original[original_site_match]
    )

    species_values[
      is.na(species_values) | trimws(species_values) == ""
    ] <- "Unassigned_species"
    site_values[
      is.na(site_values) | trimws(site_values) == ""
    ] <- "Unassigned_site"

    normalise_kinship_matrix <- function(kin, sample_names) {
      kin <- as.matrix(kin)

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

      # Account for minor numerical differences between matrix triangles.
      kin <- pmax(kin, t(kin))
      dimnames(kin) <- list(sample_names, sample_names)
      kin
    }

    sample_names <- dms_obj$sample_names

    # Scope definitions:
    #   site:
    #     1. estimate kinship once within each species using all its sites;
    #     2. mask pairs from different sites before clone-group construction.
    #   species:
    #     estimate and use all within-species relationships across sites.
    if (clone_scope == "site") {
      if (!exists("individual_kinship_by_pop_sp", mode = "function")) {
        stop(
          paste0(
            "clone_scope = 'site' requires the revised ",
            "individual_kinship_by_pop_sp() function."
          ),
          call. = FALSE
        )
      }

      kin_result <- individual_kinship_by_pop_sp(
        dart_data = dms_obj,
        basedir = RandRbase,
        species = species,
        dataset = dataset,
        pop = site_values,
        sp = species_values,
        maf = maf_val,
        mis = locus_miss,
        as_bigmat = TRUE,
        return_reference = TRUE
      )

      kin_reference <- normalise_kinship_matrix(
        kin_result$species_reference,
        sample_names
      )
      kin_decision <- normalise_kinship_matrix(
        kin_result$species_site,
        sample_names
      )

      comparison_mask <- kin_result$comparison_mask[
        sample_names,
        sample_names,
        drop = FALSE
      ]

      # A critical validation: every eligible within-site value must be
      # identical to its species-wide reference value.
      within_site_difference <- kin_decision[comparison_mask] -
        kin_reference[comparison_mask]

      if (length(within_site_difference) > 0L &&
          any(
            abs(within_site_difference) > sqrt(.Machine$double.eps),
            na.rm = TRUE
          )) {
        stop(
          paste0(
            "The site-masked matrix differs from the species-wide reference ",
            "for at least one eligible within-site pair."
          ),
          call. = FALSE
        )
      }
    } else {
      if (!exists("individual_kinship_by_pop", mode = "function")) {
        stop(
          paste0(
            "clone_scope = 'species' requires ",
            "individual_kinship_by_pop()."
          ),
          call. = FALSE
        )
      }

      kin_reference <- individual_kinship_by_pop(
        dart_data = dms_obj,
        basedir = RandRbase,
        species = species,
        dataset = dataset,
        pop = species_values,
        maf = maf_val,
        mis = locus_miss,
        as_bigmat = TRUE
      )

      kin_reference <- normalise_kinship_matrix(
        kin_reference,
        sample_names
      )
      kin_decision <- kin_reference

      comparison_mask <- outer(
        species_values,
        species_values,
        FUN = "=="
      )
      dimnames(comparison_mask) <- list(sample_names, sample_names)
    }

    list(
      decision = kin_decision,
      species_reference = kin_reference,
      comparison_mask = comparison_mask,
      species = species_values,
      site = site_values
    )
  }

  # ============================================================
  # Kinship heatmap helpers
  # ============================================================
  complete_annotation_palette <- function(values) {
    levels_needed <- unique(as.character(values[!is.na(values)]))
    levels_needed <- levels_needed[trimws(levels_needed) != ""]

    if (length(levels_needed) == 0L) {
      return(character(0))
    }

    palette_values <- colorspace::rainbow_hcl(length(levels_needed))
    names(palette_values) <- levels_needed
    palette_values
  }

  write_loop_kinship_heatmap <- function(
      kinship_matrix,
      comparison_mask,
      species_values,
      site_values,
      scope_name,
      title_text
  ) {
    sample_names <- rownames(kinship_matrix)

    if (is.null(sample_names)) {
      stop("The kinship matrix must have sample names as row names.")
    }

    plot_meta <- data.frame(
      sample = sample_names,
      Site = as.character(site_values),
      Species = as.character(species_values),
      stringsAsFactors = FALSE
    )

    if (nrow(plot_meta) != nrow(kinship_matrix)) {
      stop("Heatmap metadata and kinship matrix have different sample counts.")
    }

    sample_order <- order(
      plot_meta$Species,
      plot_meta$Site,
      plot_meta$sample,
      na.last = TRUE
    )

    ordered_samples <- plot_meta$sample[sample_order]
    plot_meta <- plot_meta[sample_order, , drop = FALSE]

    plot_matrix <- kinship_matrix[
      ordered_samples,
      ordered_samples,
      drop = FALSE
    ]

    # Grey cells are pairs excluded from clone comparison under this scope.
    # White cells are eligible comparisons with kinship near zero.
    if (!is.null(comparison_mask)) {
      ordered_mask <- comparison_mask[
        ordered_samples,
        ordered_samples,
        drop = FALSE
      ]
      plot_matrix[!ordered_mask] <- NA_real_
    }

    site_annotation <- plot_meta$Site
    site_annotation[site_annotation == "no_geo_data"] <- NA_character_
    species_annotation <- plot_meta$Species

    site_palette <- complete_annotation_palette(site_annotation)
    species_palette <- complete_annotation_palette(species_annotation)

    bottom_annotation <- ComplexHeatmap::HeatmapAnnotation(
      Site = site_annotation,
      Species = species_annotation,
      col = list(
        Site = site_palette,
        Species = species_palette
      ),
      na_col = "white",
      annotation_name_gp = grid::gpar(fontsize = 0),
      annotation_legend_param = list(
        Site = list(
          labels_gp = grid::gpar(fontsize = 6),
          title_gp = grid::gpar(fontsize = 8),
          title = site_col_name
        ),
        Species = list(
          labels_gp = grid::gpar(fontsize = 6),
          title_gp = grid::gpar(fontsize = 8),
          title = species_col_name
        )
      )
    )

    kinship_colour_function <- circlize::colorRamp2(
      c(0, 0.2, 0.4),
      c("white", "red", "black")
    )

    n_samples <- nrow(plot_matrix)
    show_sample_names <- n_samples <= kinship_show_names_max_n

    sample_name_fontsize <- if (n_samples <= 75) {
      6
    } else if (n_samples <= 150) {
      4
    } else {
      2.5
    }

    heatmap_object <- ComplexHeatmap::Heatmap(
      plot_matrix,
      col = kinship_colour_function,
      na_col = "grey92",
      bottom_annotation = bottom_annotation,
      name = "PLINK kinship",
      column_title = title_text,

      # Rasterising the matrix body is substantially faster for dense
      # kinship matrices and avoids drawing every cell as a vector object.
      use_raster = TRUE,
      raster_device = "png",
      raster_quality = 1,
      raster_by_magick = FALSE,

      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = show_sample_names,
      show_column_names = show_sample_names,
      row_names_gp = grid::gpar(fontsize = sample_name_fontsize),
      column_names_gp = grid::gpar(fontsize = sample_name_fontsize),
      row_names_max_width = grid::unit(15, "cm"),
      border_gp = grid::gpar(col = "black", lty = 1)
    )

    heatmap_file <- file.path(
      kinship_plots_dir,
      paste0("PLINK_kin_heatmap_", scope_name, ".png")
    )

    # Fixed pixel dimensions and low resolution keep diagnostic plots fast
    # and prevent very large files for datasets with many samples.
    grDevices::png(
      filename = heatmap_file,
      width = kinship_png_width_px,
      height = kinship_png_height_px,
      units = "px",
      res = kinship_png_res,
      bg = "white"
    )

    tryCatch(
      {
        ComplexHeatmap::draw(
          heatmap_object,
          merge_legend = TRUE
        )
      },
      finally = {
        grDevices::dev.off()
      }
    )

    message("Saved kinship heatmap: ", heatmap_file)
    heatmap_file
  }

  # ============================================================
  # Initial/final-only kinship output helpers
  # ============================================================
  get_locus_signature <- function(dms_obj) {
    if (!is.null(dms_obj$locus_names)) {
      return(as.character(dms_obj$locus_names))
    }

    locus_names <- colnames(dms_obj$gt)
    if (is.null(locus_names)) {
      return(as.character(seq_len(ncol(dms_obj$gt))))
    }

    as.character(locus_names)
  }

  save_initial_kinship_outputs <- function(kin_result) {
    if (isTRUE(initial_kinship_saved)) {
      return(invisible(NULL))
    }

    if (isTRUE(save_kinship_matrices)) {
      saveRDS(
        kin_result$decision,
        file = initial_decision_file,
        compress = compress_kinship_matrices
      )
      kinship_files <<- c(kinship_files, initial_decision_file)
      message("Saved initial kinship matrix: ", initial_decision_file)
    }

    if (isTRUE(write_kinship_heatmaps)) {
      initial_title <- if (clone_scope == "site") {
        paste0(
          "Species-wide kinship reference; displayed only within ",
          "species × site (initial)"
        )
      } else {
        "Kinship within species across all sites (initial)"
      }

      initial_plot <- write_loop_kinship_heatmap(
        kinship_matrix = kin_result$decision,
        comparison_mask = kin_result$comparison_mask,
        species_values = kin_result$species,
        site_values = kin_result$site,
        scope_name = paste0(clone_scope, "_initial"),
        title_text = initial_title
      )

      kinship_plot_files <<- c(
        kinship_plot_files,
        initial_plot
      )
    }

    initial_kinship_saved <<- TRUE
    invisible(NULL)
  }

  update_final_kinship_candidate <- function(
      kin_result,
      dms_obj
  ) {
    # A zero-removal clone check is a candidate final state. If the next round
    # has the same samples and loci, the resulting kinship matrix is identical,
    # so the file does not need to be written again.
    current_samples <- as.character(dms_obj$sample_names)
    current_loci <- get_locus_signature(dms_obj)

    same_as_saved_candidate <- isTRUE(
      final_kinship_candidate_ready
    ) &&
      identical(
        current_samples,
        final_candidate_samples
      ) &&
      identical(
        current_loci,
        final_candidate_loci
      )

    if (same_as_saved_candidate) {
      message(
        "Skipped duplicate final-candidate kinship output; ",
        "samples and loci were unchanged."
      )
      return(invisible(NULL))
    }

    if (isTRUE(save_kinship_matrices) ||
        isTRUE(write_kinship_heatmaps)) {
      saveRDS(
        kin_result$decision,
        file = final_candidate_cache_file,
        compress = if (isTRUE(save_kinship_matrices)) {
          compress_kinship_matrices
        } else {
          FALSE
        }
      )
    }

    if (isTRUE(save_kinship_matrices)) {
      kinship_files <<- c(
        kinship_files,
        final_decision_file
      )
      message(
        "Updated final-candidate kinship matrix: ",
        final_decision_file
      )

      # For site-level clone removal, retain one final unmasked species-wide
      # reference matrix. No intermediate reference matrices are written.
      if (clone_scope == "site" &&
          isTRUE(save_reference_kinship_matrices)) {
        saveRDS(
          kin_result$species_reference,
          file = final_reference_file,
          compress = compress_kinship_matrices
        )
        kinship_files <<- c(
          kinship_files,
          final_reference_file
        )
        message(
          "Updated final species-reference kinship matrix: ",
          final_reference_file
        )
      }
    }

    final_candidate_samples <<- current_samples
    final_candidate_loci <<- current_loci
    final_kinship_candidate_ready <<- TRUE

    invisible(NULL)
  }

  finalise_kinship_outputs <- function(final_dms) {
    outputs_requested <- isTRUE(save_kinship_matrices) ||
      isTRUE(write_kinship_heatmaps)

    if (!outputs_requested ||
        length(final_dms$sample_names) <= 1L) {
      return(invisible(NULL))
    }

    current_samples <- as.character(final_dms$sample_names)
    current_loci <- get_locus_signature(final_dms)

    candidate_is_current <- isTRUE(
      final_kinship_candidate_ready
    ) &&
      file.exists(final_candidate_cache_file) &&
      identical(
        current_samples,
        final_candidate_samples
      ) &&
      identical(
        current_loci,
        final_candidate_loci
      )

    # This fallback should rarely be needed. It guarantees that the saved
    # final matrix represents the actual final dms object if the loop ended
    # through max_rounds or another unusual path.
    if (!candidate_is_current) {
      message(
        "Recalculating kinship once for the final dms object because no ",
        "current final candidate was available."
      )

      final_result <- calculate_clone_kinship(final_dms)

      saveRDS(
        final_result$decision,
        file = final_candidate_cache_file,
        compress = if (isTRUE(save_kinship_matrices)) {
          compress_kinship_matrices
        } else {
          FALSE
        }
      )

      if (isTRUE(save_kinship_matrices)) {
        kinship_files <<- c(
          kinship_files,
          final_decision_file
        )

        if (clone_scope == "site" &&
            isTRUE(save_reference_kinship_matrices)) {
          saveRDS(
            final_result$species_reference,
            file = final_reference_file,
            compress = compress_kinship_matrices
          )
          kinship_files <<- c(
            kinship_files,
            final_reference_file
          )
        }
      }

      final_candidate_samples <<- current_samples
      final_candidate_loci <<- current_loci
      final_kinship_candidate_ready <<- TRUE
    }

    if (isTRUE(write_kinship_heatmaps)) {
      final_matrix <- readRDS(
        final_candidate_cache_file
      )

      final_matrix <- final_matrix[
        current_samples,
        current_samples,
        drop = FALSE
      ]

      final_meta <- as.data.frame(
        final_dms$meta$analyses,
        stringsAsFactors = FALSE
      )
      if (!"sample" %in% names(final_meta)) {
        final_meta$sample <- final_dms$sample_names
      }

      final_meta <- final_meta[
        match(current_samples, final_meta$sample),
        ,
        drop = FALSE
      ]

      if (anyNA(final_meta$sample)) {
        stop(
          "Final heatmap metadata could not be aligned to the final matrix.",
          call. = FALSE
        )
      }

      final_species <- as.character(
        final_meta[[species_col_name]]
      )

      final_site_match <- match(
        current_samples,
        orig_site_lookup$sample
      )
      if (anyNA(final_site_match)) {
        stop(
          "Final heatmap samples could not be matched to original sites.",
          call. = FALSE
        )
      }

      final_sites <- as.character(
        orig_site_lookup$site_original[final_site_match]
      )

      final_species[
        is.na(final_species) | trimws(final_species) == ""
      ] <- "Unassigned_species"
      final_sites[
        is.na(final_sites) | trimws(final_sites) == ""
      ] <- "Unassigned_site"

      if (clone_scope == "site") {
        final_mask <- outer(
          final_species,
          final_species,
          FUN = "=="
        ) & outer(
          final_sites,
          final_sites,
          FUN = "=="
        )

        final_title <- paste0(
          "Species-wide kinship reference; displayed only within ",
          "species × site (final)"
        )
      } else {
        final_mask <- outer(
          final_species,
          final_species,
          FUN = "=="
        )

        final_title <- "Kinship within species across all sites (final)"
      }

      dimnames(final_mask) <- list(
        current_samples,
        current_samples
      )

      final_plot <- write_loop_kinship_heatmap(
        kinship_matrix = final_matrix,
        comparison_mask = final_mask,
        species_values = final_species,
        site_values = final_sites,
        scope_name = paste0(clone_scope, "_final"),
        title_text = final_title
      )

      kinship_plot_files <<- c(
        kinship_plot_files,
        final_plot
      )
    }

    # Delete the temporary cache when matrices were not requested.
    if (!isTRUE(save_kinship_matrices) &&
        file.exists(final_candidate_cache_file)) {
      unlink(final_candidate_cache_file)
    }

    invisible(NULL)
  }

  # ============================================================
  # Identify clone components and keep the least-missing sample
  # ============================================================
  filter_clones_once <- function(dms_obj, round_number, stage = "main") {
    sample_names <- dms_obj$sample_names
    n_samples <- length(sample_names)

    if (n_samples <= 1L) {
      return(list(
        dms = dms_obj,
        removed = 0L,
        removed_samples = character(),
        assignments = data.frame()
      ))
    }

    kin_result <- calculate_clone_kinship(dms_obj)
    kin <- kin_result$decision
    kin_reference <- kin_result$species_reference

    # Retain only the first clone-decision matrix and heatmap here.
    # The final matrix is captured after a zero-removal clone check and its
    # heatmap is written once after the loop has converged.
    save_initial_kinship_outputs(kin_result)

    # Build an edge list directly from clone-like pairs. This avoids creating
    # an additional dense adjacency matrix and is much lighter for large data.
    clone_pairs <- which(
      upper.tri(kin) & kin >= clonal_threshold,
      arr.ind = TRUE
    )

    if (nrow(clone_pairs) == 0L) {
      clone_membership <- seq_len(n_samples)
      names(clone_membership) <- sample_names
    } else {
      edge_matrix <- cbind(
        sample_names[clone_pairs[, 1]],
        sample_names[clone_pairs[, 2]]
      )

      network <- igraph::graph_from_edgelist(
        edge_matrix,
        directed = FALSE
      )

      missing_vertices <- setdiff(sample_names, igraph::V(network)$name)
      if (length(missing_vertices) > 0L) {
        network <- igraph::add_vertices(
          network,
          nv = length(missing_vertices),
          name = missing_vertices
        )
      }

      clone_membership <- igraph::components(network)$membership
      clone_membership <- clone_membership[sample_names]
    }

    clones <- data.frame(
      sample = sample_names,
      genet = as.integer(unname(clone_membership[sample_names])),
      stringsAsFactors = FALSE
    )

    meta <- as.data.frame(dms_obj$meta$analyses, stringsAsFactors = FALSE)
    if (!"sample" %in% names(meta)) meta$sample <- sample_names
    meta <- meta[match(clones$sample, meta$sample), , drop = FALSE]

    missing_counts <- rowSums(is.na(dms_obj$gt))
    names(missing_counts) <- sample_names

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
    if (length(metadata_columns) > 0L) {
      clones_out <- cbind(
        clones_out,
        meta[, metadata_columns, drop = FALSE]
      )
    }

    group_sizes <- table(clones_out$genet)
    clones_out$group_size <- as.integer(
      group_sizes[as.character(clones_out$genet)]
    )
    clones_out$clone_detected <- clones_out$group_size > 1L
    clones_out$retained <- TRUE

    clone_genets <- as.integer(names(group_sizes[group_sizes > 1L]))

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
      site_match <- match(
        clones_out$sample,
        orig_site_lookup$sample
      )
      paste(
        as.character(clones_out[[species_col_name]]),
        as.character(orig_site_lookup$site_original[site_match]),
        sep = "__"
      )
    } else {
      as.character(clones_out[[species_col_name]])
    }

    clones_out$removal_reason <- ifelse(
      clones_out$clone_detected & !clones_out$retained,
      paste0("clone_removed_", clone_scope),
      ifelse(
        clones_out$clone_detected,
        "clone_representative_retained",
        "not_a_clone"
      )
    )

    clones_out <- clones_out[
      order(
        clones_out$genet,
        !clones_out$retained,
        clones_out$n_missing_loci,
        clones_out$sample
      ),
      ,
      drop = FALSE
    ]

    samples_to_keep <- clones_out$sample[clones_out$retained]
    removed_samples <- clones_out$sample[!clones_out$retained]

    if (isTRUE(write_clone_tables)) {
      table_to_write <- clones_out
      if (!isTRUE(clone_table_include_singletons)) {
        table_to_write <- table_to_write[
          table_to_write$clone_detected,
          ,
          drop = FALSE
        ]
      }

      stage_suffix <- if (identical(stage, "main")) {
        ""
      } else {
        paste0("_", stage)
      }

      output_file <- file.path(
        tables_dir,
        paste0(
          "PLINK_clones_",
          clone_scope,
          "_round",
          round_number,
          stage_suffix,
          ".xlsx"
        )
      )

      openxlsx::write.xlsx(
        table_to_write,
        file = output_file,
        asTable = FALSE,
        overwrite = TRUE
      )
    }

    # Filter the current dms object directly so no removed sample can be
    # accidentally reintroduced from an earlier intermediate object.
    filtered_obj <- remove.by.list(dms_obj, samples_to_keep)
    removed_n <- n_samples - length(filtered_obj$sample_names)

    if (removed_n < 0L) {
      stop(
        paste0(
          "Clone filtering increased the number of samples, indicating ",
          "unexpected remove.by.list behaviour."
        ),
        call. = FALSE
      )
    }

    if (removed_n != length(removed_samples)) {
      warning(
        "The number removed by remove.by.list (",
        removed_n,
        ") differs from the clone-decision table (",
        length(removed_samples),
        ")."
      )
    }

    # A zero-removal check is a potential final state. Save it under stable
    # final filenames; later equivalent checks are skipped, and changed
    # candidates overwrite the same files rather than creating duplicates.
    if (removed_n == 0L) {
      update_final_kinship_candidate(
        kin_result = kin_result,
        dms_obj = dms_obj
      )
    }

    # Do not retain the dense kinship matrix in memory after this check.
    # Initial/final output helpers have already cached any required matrix.
    list(
      dms = filtered_obj,
      removed = as.integer(removed_n),
      removed_samples = removed_samples,
      assignments = clones_out
    )
  }

  # ============================================================
  # Maximin downsampling, vectorised over the distance matrix
  # ============================================================
  downsample_max_diversity <- function(geno_mat, sample_ids, n_keep) {
    n_samples <- length(sample_ids)
    if (n_samples <= n_keep) return(sample_ids)

    d <- as.matrix(stats::dist(geno_mat))
    diag(d) <- NA_real_

    avg_dist <- rowMeans(d, na.rm = TRUE)
    avg_dist[!is.finite(avg_dist)] <- -Inf

    first_idx <- which.max(avg_dist)
    selected_idx <- first_idx

    min_distance_to_selected <- d[, first_idx]
    min_distance_to_selected[first_idx] <- -Inf
    min_distance_to_selected[
      !is.finite(min_distance_to_selected)
    ] <- -Inf

    while (length(selected_idx) < n_keep) {
      next_idx <- which.max(min_distance_to_selected)

      if (length(next_idx) == 0L ||
          !is.finite(min_distance_to_selected[next_idx])) {
        # Deterministic fallback when all remaining distances are unavailable.
        next_idx <- setdiff(seq_len(n_samples), selected_idx)[1]
      }

      selected_idx <- c(selected_idx, next_idx)

      candidate_distances <- d[, next_idx]
      candidate_distances[!is.finite(candidate_distances)] <- Inf

      current_distances <- min_distance_to_selected
      current_distances[!is.finite(current_distances)] <- Inf

      min_distance_to_selected <- pmin(
        current_distances,
        candidate_distances
      )
      min_distance_to_selected[selected_idx] <- -Inf
    }

    sample_ids[selected_idx]
  }

  # ============================================================
  # Main iterative cleaning loop
  # ============================================================
  while (TRUE) {
    round <- round + 1L

    if (round > max_rounds) {
      message("Reached max_rounds cap (", max_rounds, "). Stopping loop.")
      break
    }

    message("=== ROUND ", round, " ===")
    removed_this_round_total <- 0L

    # ----------------------------------------------------------
    # 1) Remove high-missingness samples
    # ----------------------------------------------------------
    sample_missingness <- rowMeans(is.na(dms$gt))
    names(sample_missingness) <- dms$sample_names

    samples_high_missing <- names(sample_missingness)[
      sample_missingness > sample_miss
    ]

    removed_missing <- length(samples_high_missing)

    if (removed_missing > 0L) {
      write.table(
        data.frame(
          sample = samples_high_missing,
          sample_miss = unname(sample_missingness[samples_high_missing]),
          stringsAsFactors = FALSE
        ),
        file = file.path(
          tables_dir,
          paste0("high_missing_samples_removed_round", round, ".tsv")
        ),
        sep = "\t",
        row.names = FALSE
      )

      dms <- remove.by.missingness(dms, sample_miss)
      removed_this_round_total <- removed_this_round_total + removed_missing

      message(
        "Removed ",
        removed_missing,
        " samples for high missingness (round ",
        round,
        ")."
      )
    } else {
      message("No high-missingness samples found (round ", round, ").")
    }

    add_log(round, "missingness", removed_missing)

    # Reanalyse only when the sample set changed, or when the previous round
    # ended with a clone removal after its final reanalysis.
    if (removed_missing > 0L || isTRUE(needs_reanalysis)) {
      re_res <- reanalyse(d1, dms, paste0(round, "a"))
      d1 <- re_res$d1
      dms <- re_res$dms
      treatment <- dms$treatment
      needs_reanalysis <- FALSE
    } else {
      message("Skipped unchanged missingness reanalysis.")
    }

    # ----------------------------------------------------------
    # 2) Main clone detection and removal
    # ----------------------------------------------------------
    dms <- disambiguate_site_names(dms)

    clone_result <- filter_clones_once(
      dms_obj = dms,
      round_number = round,
      stage = "main"
    )

    dms <- clone_result$dms
    removed_main_clones <- clone_result$removed
    removed_this_round_total <-
      removed_this_round_total + removed_main_clones

    add_log(
      round,
      paste0("clones_", clone_scope),
      removed_main_clones
    )

    if (removed_main_clones > 0L) {
      message(
        "Removed ",
        removed_main_clones,
        " clone samples using the '",
        clone_scope,
        "' scope (round ",
        round,
        ")."
      )
    } else {
      message(
        "No clone samples detected using the '",
        clone_scope,
        "' scope (round ",
        round,
        ")."
      )
    }

    # A clone removal can alter locus missingness and the one-SNP-per-locus
    # data. Match the original workflow by rebuilding before later filters.
    reanalysed_after_main_clone <- FALSE

    if (removed_main_clones > 0L) {
      re_res2 <- reanalyse(d1, dms, paste0(round, "b"))
      d1 <- re_res2$d1
      dms <- re_res2$dms
      treatment <- dms$treatment
      reanalysed_after_main_clone <- TRUE
    } else {
      message("Skipped unchanged post-clone reanalysis.")
    }

    dms <- disambiguate_site_names(dms)

    # ----------------------------------------------------------
    # 3) Remove small populations
    # ----------------------------------------------------------
    removed_small_pop_samples <- 0L

    if (remove_pops_flag) {
      site_values <- as.character(
        dms$meta$analyses[, site_col_name]
      )
      site_counts <- table(site_values)
      small_sites <- names(
        site_counts[site_counts < samples_per_pop_remove]
      )

      samples_to_keep <- dms$sample_names[
        !(site_values %in% small_sites)
      ]
      removed_small_pop_samples <-
        length(dms$sample_names) - length(samples_to_keep)

      if (removed_small_pop_samples > 0L) {
        dms <- remove.by.list(dms, samples_to_keep)
      }

      message(
        "Removed ",
        length(small_sites),
        " sites containing ",
        removed_small_pop_samples,
        " samples because site size was < ",
        samples_per_pop_remove,
        " (round ",
        round,
        ")."
      )
    }

    add_log(round, "remove_small_pops", removed_small_pop_samples)

    # ----------------------------------------------------------
    # 4) Downsample large populations
    # ----------------------------------------------------------
    removed_downsample <- 0L

    if (downsample_flag) {
      geno <- dms$gt
      site_values <- as.character(
        dms$meta$analyses[, site_col_name]
      )
      sites <- unique(site_values)
      samples_to_keep <- character()

      for (site_name in sites) {
        site_idx <- which(site_values == site_name)
        site_samples <- dms$sample_names[site_idx]

        if (length(site_samples) > samples_per_pop_downsample) {
          geno_sub <- geno[site_idx, , drop = FALSE]

          keep <- downsample_max_diversity(
            geno_mat = geno_sub,
            sample_ids = site_samples,
            n_keep = samples_per_pop_downsample
          )
          samples_to_keep <- c(samples_to_keep, keep)
        } else {
          samples_to_keep <- c(samples_to_keep, site_samples)
        }
      }

      samples_to_keep <- unique(samples_to_keep)
      removed_downsample <-
        length(dms$sample_names) - length(samples_to_keep)

      if (removed_downsample > 0L) {
        dms <- remove.by.list(dms, samples_to_keep)
      }

      message(
        "Downsampling removed ",
        removed_downsample,
        " samples and retained at most ",
        samples_per_pop_downsample,
        " per site (round ",
        round,
        ")."
      )
    }

    add_log(round, "downsample", removed_downsample)

    removed_pop_total <-
      removed_small_pop_samples + removed_downsample
    removed_this_round_total <-
      removed_this_round_total + removed_pop_total

    # Rebuild only if population filtering or downsampling changed samples.
    if (removed_pop_total > 0L) {
      re_res3 <- reanalyse(d1, dms, paste0(round, "c"))
      d1 <- re_res3$d1
      dms <- re_res3$dms
      treatment <- dms$treatment
      reanalysed_after_main_clone <- TRUE
    } else {
      message("Skipped unchanged population/downsampling reanalysis.")
    }

    # ----------------------------------------------------------
    # 5) Final clone check
    # ----------------------------------------------------------
    # This check is necessary only when loci were recalculated after the main
    # clone check. If no reanalysis occurred, it would repeat the same matrix.
    removed_final_clones <- 0L

    if (isTRUE(reanalysed_after_main_clone)) {
      dms <- disambiguate_site_names(dms)

      clone_final_result <- filter_clones_once(
        dms_obj = dms,
        round_number = round,
        stage = "finalcheck"
      )

      dms <- clone_final_result$dms
      removed_final_clones <- clone_final_result$removed
      removed_this_round_total <-
        removed_this_round_total + removed_final_clones

      if (removed_final_clones > 0L) {
        needs_reanalysis <- TRUE
        message(
          "Final clone check removed ",
          removed_final_clones,
          " samples using the '",
          clone_scope,
          "' scope (round ",
          round,
          ")."
        )
      } else {
        message(
          "No clones detected in the final '",
          clone_scope,
          "' check (round ",
          round,
          ")."
        )
      }
    } else {
      message(
        "Skipped final clone check because loci were unchanged after the ",
        "main clone check."
      )
    }

    add_log(
      round,
      paste0("clones_", clone_scope, "_finalcheck"),
      removed_final_clones
    )

    # ----------------------------------------------------------
    # Round summary and stopping rule
    # ----------------------------------------------------------
    new_n <- length(dms$sample_names)

    if (isTRUE(write_round_summaries)) {
      round_summary <- data.frame(
        round = round,
        samples_after = new_n,
        removed_this_round = removed_this_round_total,
        stringsAsFactors = FALSE
      )

      write.table(
        round_summary,
        file = file.path(
          tables_dir,
          paste0("round_summary_", round, ".tsv")
        ),
        sep = "\t",
        row.names = FALSE
      )
    }

    message(
      "Round ",
      round,
      " end: samples = ",
      new_n,
      " (previous ",
      prev_n,
      ")."
    )

    if (new_n == prev_n) {
      message("No change in sample count; stopping.")
      break
    }

    prev_n <- new_n
  }

  # ============================================================
  # Restore original site names
  # ============================================================
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

  # Write one final heatmap and ensure the final matrix files correspond
  # exactly to the converged dms object.
  finalise_kinship_outputs(dms)

  log_df <- if (length(log_rows) > 0L) {
    do.call(rbind, log_rows)
  } else {
    data.frame(
      round = integer(),
      step = character(),
      removed = integer(),
      stringsAsFactors = FALSE
    )
  }

  list(
    dms = dms,
    d1 = d1,
    rounds = round,
    log = log_df,
    treatment = treatment,
    clone_scope = clone_scope,
    kinship_files = unique(kinship_files),
    kinship_plot_files = unique(kinship_plot_files),
    optimisation_settings = list(
      write_intermediate_qc = write_intermediate_qc,
      write_clone_tables = write_clone_tables,
      clone_table_include_singletons = clone_table_include_singletons,
      write_round_summaries = write_round_summaries,
      save_kinship_matrices = save_kinship_matrices,
      save_reference_kinship_matrices = save_reference_kinship_matrices,
      write_kinship_heatmaps = write_kinship_heatmaps,
      kinship_show_names_max_n = kinship_show_names_max_n,
      kinship_png_width_px = kinship_png_width_px,
      kinship_png_height_px = kinship_png_height_px,
      kinship_png_res = kinship_png_res,
      matrix_format = "rds",
      matrix_compression = compress_kinship_matrices,
      matrix_output_mode = "initial_final_only",
      heatmap_format = "png",
      heatmap_rasterised = TRUE,
      heatmap_output_mode = "initial_final_only",
      final_species_reference_saved = (
        clone_scope == "site" &&
          isTRUE(save_kinship_matrices) &&
          isTRUE(save_reference_kinship_matrices)
      ),
      kinship_directory = if (isTRUE(save_kinship_matrices)) kinship_dir else NULL,
      kinship_plot_directory = if (isTRUE(write_kinship_heatmaps)) kinship_plots_dir else NULL
    )
  )
}
