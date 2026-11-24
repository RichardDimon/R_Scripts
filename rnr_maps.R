###############################################################################
# rnr_maps.R
# Sampling-site maps for the RnR dms analysis workflow
# ("1. RnR dms analysis workflow.R").
#
# Source it from the main workflow, in the place where the MAPS section was:
#   source_url_once(url = "https://github.com/RichardDimon/R_Scripts/blob/main/rnr_maps.R?raw=TRUE",
#                   script_name = "rnr_maps", force_reload = reload_github_functions)
#   make_rnr_base_map()                       # always: base_map/divxlims are used by the later scatterpie maps
#   run_rnr_maps(maps = c("site", "species")) # the maps to draw
#
# Sourcing only defines settings and functions, so it takes about a second.
# Maps are drawn only by run_rnr_maps(), which runs only the preparation steps
# that the requested maps need.
#
# Requires (made earlier in the main workflow): filtered_site_summary,
# output_dir, site_col_name, species_col_name.
# Optional: site_colours, sp_colours, sp_shapes, fullspecies.
#
# Maps (use the name or the number in run_rnr_maps(maps = ...)):
#   0  "context"            Map_00_Australia_context
#   1  "site_local"         Map_01_Sites_local_labels
#   2  "site_ordered"       Map_02_Sites_ordered_labels
#   3  "species_local"      Map_03_Species_local_labels
#   4  "species_ordered"    Map_04_Species_ordered_labels
#   5  "species_no_labels"  Map_05_Species_no_labels
#   6  "sampling_effort"    Map_06_Sampling_effort
#   7  "species_facets"     Map_07_Species_facets (only when >1 species)
#   8  "longitude_profile"  Map_08_Longitude_profile
#   9  "latitude_profile"   Map_09_Latitude_profile
#   10 "species_points"     Map_10_All_species_points_with_inset
#   11 "site"               Map_11_All_site_labels_with_inset              (default)
#   12 "species"            Map_12_All_site_labels_by_species_with_inset   (default)
#   13 "site_numbers"       Map_13_All_site_numbers_with_inset
#   14 "species_numbers"    Map_14_All_site_numbers_by_species_with_inset
#   "all" = all of the above; "none" (or FALSE / character(0)) = no maps.
#
# Any setting in section 1 (map_dpi, map_save_pdf, map_inset_*, ...) can be
# changed by assigning it in the main workflow before calling run_rnr_maps().
###############################################################################

# Environment the script was sourced into (.GlobalEnv from source_url_once()).
rnr_maps_env <- environment()


#### 1. Default map settings ####
# A setting that already exists in the session (for example set in the main
# workflow) is kept; only missing settings are filled with these defaults.
local({
  # --- User-adjustable map settings
  map_save_pdf <- FALSE # TRUE also saves PDFs
  map_dpi <- 600
  map_repel_label_cap <- 120
  map_facet_label_cap <- 14
  map_facet_max_columns <- 5L # automatic facet-layout ceiling
  map_label_cluster_distance <- 0.14 # larger groups more nearby labels
  map_label_radial_spread <- 1.65 # larger uses more map space
  map_min_detail_aspect <- 0.72
  map_min_detail_width_cm <- 32
  map_seed <- 417
  # ADDITIONAL SPECIES-COLOURED AND SITE-LABELLED MAPS WITH AUTOMATIC INSETS
  #

  #  User-adjustable inset settings

  map_inset_radius_km <- NULL
  map_inset_min_sites <- 5L
  map_inset_max_clusters <- 3L
  map_inset_halo_steps <- 1L
  # An inset is only required when at least this many sites participate in close
  # screen-space links. This prevents sparse datasets from receiving an inset
  # merely because four sites happen to fall within the geographic search radius.
  map_inset_visual_trigger_prop <- 0.030
  map_inset_min_close_pairs <- 5L
  # Additional inset panels have a stricter threshold than the primary inset.
  # This prevents a weak secondary grouping from creating another large panel.
  map_inset_secondary_min_sites <- 6L
  map_inset_secondary_min_close_pairs <- 8L
  map_inset_secondary_density_ratio <- 0.55
  map_inset_cluster_separation_factor <- 1.75
  map_inset_padding <- 0.10
  map_inset_aspect_fit_strength <- 1.00
  map_inset_width <- 0.30
  map_inset_height <- 0.34
  # Every inset is resized independently from the geographic aspect of its own
  # cluster. These bounds prevent either a very narrow portrait or a very wide
  # landscape inset.
  map_inset_min_frame_aspect <- 2 / 3
  map_inset_max_frame_aspect <- 3 / 2
  map_inset_min_fraction <- 0.18
  map_inset_max_fraction <- 0.46
  # Inset frames can grow when a large proportion of the map is available ocean
  # or point-free space. Growth is capped so the inset never dominates the map.
  map_inset_max_adaptive_scale <- 1.40
  map_inset_blank_area_reference <- 0.25
  # If a full inset does not initially fit over water, inspect the space just
  # outside each map edge and expand toward a strongly ocean-dominated side before
  # considering an inland fallback.
  map_inset_ocean_detection_threshold <- 0.65
  map_inset_ocean_expansion_step <- 0.10
  map_inset_ocean_expansion_attempts <- 7L
  map_inset_outer_margin <- 0.012
  map_inset_extra_side_prop <- 0.10
  map_inset_extra_vertical_prop <- 0.025
  map_outer_padding_fraction <- 0.035
  map_outer_minimum_padding <- 0.045
  # Candidate corners and edge-centred positions are ranked jointly. Covering a
  # mapped site is penalised most strongly, followed by connector crossings and
  # land coverage; distance from the dense cluster is retained only as a mild
  # aesthetic tie-breaker.
  map_inset_land_grid_size <- 14L
  map_inset_site_penalty <- 100
  map_inset_connector_penalty <- 12
  map_inset_land_penalty <- 10
  map_inset_distance_penalty <- 0.35
  map_inset_connector_buffer <- 0.022
  map_inset_overlap_penalty <- 1000
  # Require a visible gap between neighbouring inset frames. This is especially
  # useful when two insets occupy the same side of a map: the candidate selector
  # will favour a more distant edge/corner position instead of stacking the
  # frames tightly together.
  map_inset_min_spacing <- 0.065
  # A side is treated as ocean-facing when the proposed inset footprints contain
  # at least this proportion of water. Multi-inset maps preferentially use one
  # such side and distribute the frames evenly along it.
  map_inset_ocean_stack_threshold <- 0.55
  # All outputs produced by the MAPS section are kept between 2:3 portrait and
  # 3:2 landscape. The shorter device dimension is enlarged when necessary.
  map_output_min_aspect <- 2 / 3
  map_output_max_aspect <- 3 / 2
  map_labelled_point_size <- 4.2
  map_inset_point_scale <- 1
  map_inset_label_text_size <- 3.5
  # Numbered-site alternatives use larger text because the labels are short.
  # These values control Maps 13-14 independently of the full-name maps.
  map_number_label_text_size <- 7.4
  map_number_key_text_size <- 6.2
  map_number_key_point_size <- 6.1
  # Long key entries wrap onto a second line instead of extending into the next
  # column. The row spacing accommodates two lines at this width.
  map_number_key_wrap_width <- 34L
  # Leave enough separation for the enlarged two-column key so long site names
  # cannot collide with the marker or number at the start of column two.
  map_number_key_column_spacing <- 1.16
  map_species_legend_point_size <- 5
  map_labelled_dpi <- suppressWarnings(as.integer(
    getOption("map.labelled.dpi", 800L)
  ))

  if (!is.finite(map_labelled_dpi) || map_labelled_dpi < 72L) {
    map_labelled_dpi <- 800L
  }

  for (setting_name in setdiff(ls(environment()), "setting_name")) {
    if (!exists(setting_name, envir = rnr_maps_env, inherits = FALSE)) {
      assign(setting_name, get(setting_name, inherits = FALSE), envir = rnr_maps_env)
    }
  }
})


#### 2. Map helper functions ####
# Defined when the script is sourced, so later sections (for example the
# fastDiversity maps, which reuse add_labelled_map_inset(),
# score_map_inset_candidates() and constrain_map_output_dimensions()) work
# even when no maps are drawn.

# Automatic map limits
make_auto_map_limits <- function(df, x_col = "long", y_col = "lat", pad_prop = 0.055,
                                 min_pad_x = 0.08, min_pad_y = 0.08, target_aspect = 1.30) {
  x_rng <- range(df[[x_col]], na.rm = TRUE); y_rng <- range(df[[y_col]], na.rm = TRUE)
  x_span <- diff(x_rng); y_span <- diff(y_rng)
  if (!is.finite(x_span) || x_span == 0) x_span <- min_pad_x * 2
  if (!is.finite(y_span) || y_span == 0) y_span <- min_pad_y * 2
  x_pad <- max(x_span * pad_prop, min_pad_x); y_pad <- max(y_span * pad_prop, min_pad_y)
  xlim <- c(x_rng[1] - x_pad, x_rng[2] + x_pad); ylim <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)
  x_span2 <- diff(xlim); y_span2 <- diff(ylim); current_aspect <- x_span2 / y_span2
  if (is.finite(current_aspect) && current_aspect < 1 / target_aspect) {
    extra <- (y_span2 / target_aspect - x_span2) / 2
    if (extra > 0) xlim <- xlim + c(-extra, extra)
  }
  if (is.finite(current_aspect) && current_aspect > target_aspect) {
    extra <- (x_span2 / target_aspect - y_span2) / 2
    if (extra > 0) ylim <- ylim + c(-extra, extra)
  }
  list(xlim = xlim, ylim = ylim)
}

# Single right-hand label column
make_single_column_label_positions <- function(label_df, xlim, ylim, side = "right", x_gap_prop = 0.025,
                                               label_margin_prop = NULL, top_bottom_pad_prop = 0.012,
                                               blend_actual_lat = 0.12) {
  label_df <- label_df %>% dplyr::arrange(dplyr::desc(lat))
  n <- nrow(label_df); x_span <- diff(xlim); y_span <- diff(ylim)
  y_top <- ylim[2] - y_span * top_bottom_pad_prop; y_bot <- ylim[1] + y_span * top_bottom_pad_prop
  label_df$latlabel_even <- seq(y_top, y_bot, length.out = n)
  lat_scaled <- scales::rescale(label_df$lat, to = c(y_bot, y_top))
  label_df$latlabel <- (1 - blend_actual_lat) * label_df$latlabel_even + blend_actual_lat * lat_scaled
  label_df <- label_df %>% dplyr::arrange(dplyr::desc(latlabel))
  max_label_chars <- max(nchar(as.character(label_df$site_label)), na.rm = TRUE)
  if (is.null(label_margin_prop)) label_margin_prop <- dplyr::case_when(
    max_label_chars <= 4 ~ 0.07, max_label_chars <= 10 ~ 0.12,
    max_label_chars <= 20 ~ 0.20, max_label_chars <= 35 ~ 0.30, TRUE ~ 0.44)
  if (side == "right") {
    label_df$longlabel <- max(label_df$long, na.rm = TRUE) + x_span * x_gap_prop
    xlim[2] <- max(xlim[2], label_df$longlabel[1] + x_span * label_margin_prop); hjust_val <- 0
  } else {
    label_df$longlabel <- min(label_df$long, na.rm = TRUE) - x_span * x_gap_prop
    xlim[1] <- min(xlim[1], label_df$longlabel[1] - x_span * label_margin_prop); hjust_val <- 1
  }
  list(labels = label_df, xlim = xlim, ylim = ylim, hjust = hjust_val,
       max_label_chars = max_label_chars, label_margin_prop = label_margin_prop)
}

parse_map_coordinate <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  has_south <- grepl("[Ss]", x_chr)
  has_west <- grepl("[Ww]", x_chr)
  has_north <- grepl("[Nn]", x_chr)
  has_east <- grepl("[Ee]", x_chr)
  decimal_comma <- grepl("^[[:space:]]*[+-]?[0-9]+,[0-9]+[[:space:]]*[NSEWnsew]?[[:space:]]*$", x_chr) & !grepl("\\.", x_chr)
  x_chr[decimal_comma] <- sub(",", ".", x_chr[decimal_comma], fixed = TRUE)
  numeric_text <- gsub("[NSEWnsew°º[:space:]]", "", x_chr)
  value <- suppressWarnings(as.numeric(numeric_text))
  value[has_south | has_west] <- -abs(value[has_south | has_west])
  value[has_north | has_east] <- abs(value[has_north | has_east])
  value
}

clean_map_label <- function(x) {
  x <- gsub("_+", " ", as.character(x))
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

format_sample_count <- function(x) {
  formatC(as.integer(round(x)), format = "d", big.mark = ",")
}

complete_map_palette <- function(groups, supplied = NULL, chroma = 70, luminance = 55) {
  groups <- unique(as.character(groups))
  n_groups <- length(groups)
  default_colours <- grDevices::hcl(h = seq(15, 375, length.out = n_groups + 1L)[seq_len(n_groups)], c = chroma, l = luminance)
  output <- stats::setNames(default_colours, groups)
  if (!is.null(supplied) && length(supplied) > 0) {
    supplied <- as.character(supplied)
    if (is.null(names(supplied))) {
      n_copy <- min(length(supplied), n_groups)
      output[seq_len(n_copy)] <- supplied[seq_len(n_copy)]
    }
    else {
      shared <- intersect(groups, names(supplied))
      output[shared] <- supplied[shared]
    }
  }
  valid_colour <- vapply(output, function(z) {
    !inherits(try(grDevices::col2rgb(z), silent = TRUE), "try-error")
  }, logical(1))
  output[!valid_colour | is.na(output) | !nzchar(output)] <- default_colours[!valid_colour | is.na(output) | !nzchar(output)]
  output <- scales::alpha(output, 1)
  output
}

complete_map_shapes <- function(groups, supplied = NULL) {
  groups <- unique(as.character(groups))
  default_shapes <- rep(c(21:25, 16, 17, 15, 18, 8, 3, 4, 7, 9, 10, 12), length.out = length(groups))
  output <- stats::setNames(default_shapes, groups)
  if (!is.null(supplied) && length(supplied) > 0) {
    supplied_names <- names(supplied)
    supplied <- suppressWarnings(as.numeric(supplied))
    names(supplied) <- supplied_names
    if (is.null(supplied_names)) {
      n_copy <- min(length(supplied), length(groups))
      output[seq_len(n_copy)] <- supplied[seq_len(n_copy)]
    }
    else {
      supplied_names <- names(supplied)
      shared <- intersect(groups, supplied_names)
      output[shared] <- supplied[match(shared, supplied_names)]
    }
  }
  output[!is.finite(output)] <- default_shapes[!is.finite(output)]
  output
}

make_map_limits <- function(df, pad_fraction = 0.1, minimum_span = 0.36) {
  x_range <- range(df$long, na.rm = TRUE)
  y_range <- range(df$lat, na.rm = TRUE)
  x_mid <- mean(x_range)
  y_mid <- mean(y_range)
  x_span <- max(diff(x_range), minimum_span)
  y_span <- max(diff(y_range), minimum_span)
  x_range <- x_mid + c(-0.5, 0.5) * x_span
  y_range <- y_mid + c(-0.5, 0.5) * y_span
  x_pad <- max(0.08, x_span * pad_fraction)
  y_pad <- max(0.08, y_span * pad_fraction)
  xlim <- x_range + c(-x_pad, x_pad)
  ylim <- y_range + c(-y_pad, y_pad)
  xlim <- pmax(94.5, pmin(170.5, xlim))
  ylim <- pmax(-56.5, pmin(-7.5, ylim))
  list(xlim = xlim, ylim = ylim)
}

ensure_minimum_map_aspect <- function(limits, minimum_aspect = 0.72) {
  mean_lat <- mean(limits$ylim)
  latitude_correction <- max(0.25, cos(mean_lat * pi/180))
  current_aspect <- diff(limits$xlim) * latitude_correction/diff(limits$ylim)
  if (is.finite(current_aspect) && current_aspect < minimum_aspect) {
    required_x_span <- diff(limits$ylim) * minimum_aspect/latitude_correction
    x_mid <- mean(limits$xlim)
    new_xlim <- x_mid + c(-0.5, 0.5) * required_x_span
    if (new_xlim[1] < 94.5)
      new_xlim <- new_xlim + (94.5 - new_xlim[1])
    if (new_xlim[2] > 170.5)
      new_xlim <- new_xlim - (new_xlim[2] - 170.5)
    limits$xlim <- new_xlim
  }
  limits
}

format_longitude <- function(x) {
  paste0(format(round(abs(x), 2), trim = TRUE, scientific = FALSE), "°", ifelse(x < 0, "W", "E"))
}

format_latitude <- function(x) {
  paste0(format(round(abs(x), 2), trim = TRUE, scientific = FALSE), "°", ifelse(x < 0, "S", "N"))
}

constrain_map_geographic_limits <- function(limits,
                                            minimum_aspect = 2 / 3,
                                            maximum_aspect = 3 / 2,
                                            horizontal_side = "balanced",
                                            vertical_side = "balanced") {
  xlim <- suppressWarnings(as.numeric(limits$xlim))
  ylim <- suppressWarnings(as.numeric(limits$ylim))
  if (length(xlim) != 2L || length(ylim) != 2L ||
      !all(is.finite(c(xlim, ylim))) ||
      diff(xlim) <= 0 || diff(ylim) <= 0) {
    return(limits)
  }

  longitude_scale <- max(0.20, cos(mean(ylim) * pi / 180))
  geographic_aspect <- diff(xlim) * longitude_scale / diff(ylim)

  if (geographic_aspect < minimum_aspect) {
    required_x_span <- minimum_aspect * diff(ylim) / longitude_scale
    extra_x <- max(0, required_x_span - diff(xlim))
    if (horizontal_side == "left") {
      xlim[1] <- xlim[1] - extra_x
    } else if (horizontal_side == "right") {
      xlim[2] <- xlim[2] + extra_x
    } else {
      xlim <- xlim + c(-0.5, 0.5) * extra_x
    }
  } else if (geographic_aspect > maximum_aspect) {
    required_y_span <- diff(xlim) * longitude_scale / maximum_aspect
    extra_y <- max(0, required_y_span - diff(ylim))
    if (vertical_side == "bottom") {
      ylim[1] <- ylim[1] - extra_y
    } else if (vertical_side == "top") {
      ylim[2] <- ylim[2] + extra_y
    } else {
      ylim <- ylim + c(-0.5, 0.5) * extra_y
    }
  }

  list(xlim = xlim, ylim = ylim)
}

make_map_base <- function(limits, title, subtitle = NULL, caption = NULL, ornaments = TRUE, projection_crs = 3577) {
  limits <- constrain_map_geographic_limits(limits)
  p <- ggplot2::ggplot() + ggplot2::geom_sf(data = aus_states, fill = "#F4F0E6", colour = "#6F7C82", linewidth = 0.35) + ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 5), labels = format_longitude, expand = ggplot2::expansion(mult = 0)) + ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 5), labels = format_latitude, expand = ggplot2::expansion(mult = 0)) + ggplot2::coord_sf(crs = sf::st_crs(projection_crs), default_crs = sf::st_crs(4326), datum = sf::st_crs(4326),
                                                                                                                                                                                                                                                                                                                                                                                                                        xlim = limits$xlim, ylim = limits$ylim, expand = FALSE, clip = "on") + ggplot2::labs(title = title, subtitle = subtitle, x = "Longitude", y = "Latitude", caption = caption) + ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#EAF3F5", colour = NA), panel.grid.major = ggplot2::element_line(colour = scales::alpha("white", 0.9), linewidth = 0.4), panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(fill = NA,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           colour = "#A7B1B5", linewidth = 0.45), axis.text = ggplot2::element_text(colour = "#48565C", size = 8), axis.title = ggplot2::element_text(colour = "#334148", size = 9, face = "bold"), plot.title = ggplot2::element_text(colour = "#17252B", size = 15, face = "bold"), plot.subtitle = ggplot2::element_text(colour = "#48565C", size = 9, margin = ggplot2::margin(b = 6)), plot.caption = ggplot2::element_text(colour = "#65737A", size = 7, hjust = 0), plot.background = ggplot2::element_rect(fill = "white",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   colour = NA), legend.background = ggplot2::element_rect(fill = "white", colour = NA), legend.key = ggplot2::element_rect(fill = "white", colour = NA), legend.title = ggplot2::element_text(face = "bold", size = 8), legend.text = ggplot2::element_text(size = 7), plot.margin = ggplot2::margin(8, 10, 8, 8))
  if (ornaments) {
    geographic_width <- diff(limits$xlim) * max(0.25, cos(mean(limits$ylim) * pi/180))
    scale_bar_is_stable <- max(geographic_width, diff(limits$ylim)) <= 7
    if (scale_bar_is_stable) {
      p <- p + ggspatial::annotation_scale(location = "bl", width_hint = 0.18, unit_category = "metric", style = "ticks", line_col = "#334148", text_col = "#334148", line_width = 0.55, text_cex = 0.62, pad_x = grid::unit(0.35, "cm"), pad_y = grid::unit(0.35, "cm"))
    }
    p <- p + ggspatial::annotation_north_arrow(location = "tl", which_north = "true", height = grid::unit(0.75, "cm"), width = grid::unit(0.75, "cm"), pad_x = grid::unit(0.35, "cm"), pad_y = grid::unit(0.35, "cm"), style = ggspatial::north_arrow_minimal(line_col = "#334148", fill = "#334148", text_col = "#334148", text_size = 7))
  }
  p
}

clamp_map_value <- function(x, lower, upper) {
  max(lower, min(upper, x))
}

constrain_map_output_dimensions <- function(dimensions,
                                            minimum_aspect = NULL,
                                            maximum_aspect = NULL) {
  width_cm <- suppressWarnings(as.numeric(dimensions["width"]))
  height_cm <- suppressWarnings(as.numeric(dimensions["height"]))
  if (!is.finite(width_cm) || !is.finite(height_cm) ||
      width_cm <= 0 || height_cm <= 0) {
    return(dimensions)
  }

  if (is.null(minimum_aspect)) {
    minimum_aspect <- get0(
      "map_output_min_aspect",
      envir = .GlobalEnv,
      inherits = TRUE,
      ifnotfound = 2 / 3
    )
  }
  if (is.null(maximum_aspect)) {
    maximum_aspect <- get0(
      "map_output_max_aspect",
      envir = .GlobalEnv,
      inherits = TRUE,
      ifnotfound = 3 / 2
    )
  }
  minimum_aspect <- suppressWarnings(as.numeric(minimum_aspect))[1]
  maximum_aspect <- suppressWarnings(as.numeric(maximum_aspect))[1]
  if (!is.finite(minimum_aspect) || minimum_aspect <= 0) {
    minimum_aspect <- 2 / 3
  }
  if (!is.finite(maximum_aspect) ||
      maximum_aspect < minimum_aspect) {
    maximum_aspect <- 3 / 2
  }

  output_aspect <- width_cm / height_cm
  if (output_aspect < minimum_aspect) {
    width_cm <- height_cm * minimum_aspect
  } else if (output_aspect > maximum_aspect) {
    height_cm <- width_cm / maximum_aspect
  }

  c(width = width_cm, height = height_cm)
}

map_dimensions <- function(limits, n_labels = 0, ordered_labels = FALSE) {
  mean_lat <- mean(limits$ylim)
  geographic_aspect <- diff(limits$xlim) * max(0.25, cos(mean_lat * pi/180))/diff(limits$ylim)
  geographic_aspect <- clamp_map_value(geographic_aspect, 0.18, 4.5)
  width_cm <- clamp_map_value(24 * sqrt(geographic_aspect), 18, 38)
  height_cm <- clamp_map_value(24/sqrt(geographic_aspect), 18, 44)
  if (n_labels > 35) {
    label_multiplier <- min(1.45, sqrt(n_labels/35))
    width_cm <- min(46, width_cm * sqrt(label_multiplier))
    height_cm <- min(52, height_cm * sqrt(label_multiplier))
  }
  if (ordered_labels) {
    height_cm <- max(height_cm, min(80, 11 + 0.46 * n_labels))
    width_cm <- min(52, max(width_cm, 28))
  }
  c(width = width_cm, height = height_cm)
}

select_spatial_labels <- function(df, cap) {
  if (!is.finite(cap) || nrow(df) <= cap)
    return(df)
  cap <- max(3L, as.integer(cap))
  n_each <- max(1L, floor(cap/3))
  latitude_order <- order(df$lat, df$long)
  longitude_order <- order(df$long, df$lat)
  population_order <- order(df$plot_n, decreasing = TRUE)
  spaced_indices <- function(index, n) {
    index[unique(round(seq(1, length(index), length.out = min(n, length(index)))))]
  }
  selected <- unique(c(spaced_indices(latitude_order, n_each), spaced_indices(longitude_order, n_each), head(population_order, n_each)))
  if (length(selected) < cap) {
    selected <- c(selected, head(setdiff(population_order, selected), cap - length(selected)))
  }
  df[selected[seq_len(min(cap, length(selected)))], , drop = FALSE]
}

add_local_label_nudges <- function(df, limits) {
  n_labels <- nrow(df)
  if (n_labels == 0) {
    df$label_nudge_x <- numeric()
    df$label_nudge_y <- numeric()
    return(df)
  }
  x_span <- diff(limits$xlim)
  y_span <- diff(limits$ylim)
  latitude_correction <- max(0.25, cos(mean(limits$ylim) * pi/180))
  x_scaled <- ((df$long - limits$xlim[1])/x_span) * latitude_correction
  y_scaled <- (df$lat - limits$ylim[1])/y_span
  scaled_coordinates <- cbind(x_scaled, y_scaled)
  if (n_labels == 1) {
    cluster_id <- 1L
  }
  else {
    cluster_id <- stats::cutree(stats::hclust(stats::dist(scaled_coordinates), method = "complete"), h = map_label_cluster_distance)
  }
  df$label_nudge_x <- 0
  df$label_nudge_y <- 0
  labels_per_ring <- 6L
  map_centre <- c(0.5 * latitude_correction, 0.5)
  for (cluster in unique(cluster_id)) {
    indices <- which(cluster_id == cluster)
    cluster_centre <- colMeans(scaled_coordinates[indices, , drop = FALSE])
    outward_angle <- atan2(cluster_centre[2] - map_centre[2], cluster_centre[1] - map_centre[1])
    if (!is.finite(outward_angle)) {
      outward_angle <- (cluster * 2.399963)%%(2 * pi)
    }
    if (length(indices) > 1) {
      point_angles <- atan2(scaled_coordinates[indices, 2] - cluster_centre[2], scaled_coordinates[indices, 1] - cluster_centre[1])
      label_key <- if ("display_species_site_n" %in% names(df)) {
        df$display_species_site_n[indices]
      }
      else if ("display_site_n" %in% names(df)) {
        df$display_site_n[indices]
      }
      else {
        as.character(indices)
      }
      indices <- indices[order(point_angles, label_key)]
    }
    cluster_size <- length(indices)
    for (position in seq_along(indices)) {
      ring <- (position - 1L)%/%labels_per_ring
      first_on_ring <- ring * labels_per_ring + 1L
      number_on_ring <- min(labels_per_ring, cluster_size - first_on_ring + 1L)
      slot <- position - first_on_ring
      angle <- outward_angle + 2 * pi * slot/number_on_ring + ring * pi/labels_per_ring
      radius <- min(0.4, map_label_radial_spread * (0.06 + 0.045 * ring))
      df$label_nudge_x[indices[position]] <- cos(angle) * radius * x_span
      df$label_nudge_y[indices[position]] <- sin(angle) * radius * y_span * 1.35
    }
  }
  reflect_into_interval <- function(x, lower, upper) {
    span <- upper - lower
    folded <- (x - lower)%%(2 * span)
    lower + ifelse(folded <= span, folded, 2 * span - folded)
  }
  target_x <- reflect_into_interval(df$long + df$label_nudge_x, repel_label_xlim[1], repel_label_xlim[2])
  target_y <- reflect_into_interval(df$lat + df$label_nudge_y, repel_label_ylim[1], repel_label_ylim[2])
  df$label_nudge_x <- target_x - df$long
  df$label_nudge_y <- target_y - df$lat
  df
}

make_local_label_layer <- function(data, label_field, colour_field = NULL, text_size = 3, neutral_colour = "#263238") {
  mapping <- if (is.null(colour_field)) {
    ggplot2::aes(x = .data$long, y = .data$lat, label = .data[[label_field]])
  }
  else {
    ggplot2::aes(x = .data$long, y = .data$lat, label = .data[[label_field]], colour = .data[[colour_field]])
  }
  arguments <- list(data = data, mapping = mapping, seed = map_seed, size = text_size, fontface = "bold", fill = scales::alpha("white", 0.97), label.padding = grid::unit(0.16, "lines"), box.padding = grid::unit(0.9, "lines"), point.padding = grid::unit(0.42, "lines"), min.segment.length = 0, segment.alpha = 0.88, segment.size = 0.42, force = 3.1, force_pull = 0, nudge_x = data$label_nudge_x, nudge_y = data$label_nudge_y, max.time = 40, max.iter = 150000, max.overlaps = Inf, xlim = repel_label_xlim,
                    ylim = repel_label_ylim, show.legend = FALSE)
  if (is.null(colour_field))
    arguments$colour <- neutral_colour
  do.call(ggrepel::geom_label_repel, arguments)
}

make_ordered_label_column <- function(site_df, limits, text_size) {
  output <- dplyr::arrange(site_df, dplyr::desc(.data$lat), .data$long)
  n_labels <- nrow(output)
  x_span <- diff(limits$xlim)
  y_span <- diff(limits$ylim)
  y_top <- limits$ylim[2] - 0.025 * y_span
  y_bottom <- limits$ylim[1] + 0.025 * y_span
  if (n_labels == 1) {
    y_even <- mean(c(y_top, y_bottom))
    y_geographic <- y_even
  }
  else {
    y_even <- seq(y_top, y_bottom, length.out = n_labels)
    y_geographic <- scales::rescale(output$lat, to = c(y_bottom, y_top), from = range(output$lat))
    if (!all(is.finite(y_geographic)))
      y_geographic <- y_even
  }
  output$label_lat <- 0.78 * y_even + 0.22 * y_geographic
  output$label_long <- limits$xlim[2] + 0.045 * x_span
  output$label_text_long <- output$label_long + 0.018 * x_span
  max_characters <- max(nchar(output$display_site), na.rm = TRUE)
  label_space <- max(0.35, min(6, x_span * (0.035 + 0.0065 * max_characters * text_size/3)))
  right_margin <- max(0.1, 0.05 * label_space)
  expanded_limits <- limits
  expanded_limits$xlim[2] <- max(expanded_limits$xlim[2], output$label_text_long[1] + label_space + right_margin)
  list(data = output, limits = expanded_limits)
}

add_sample_size_scale <- function(p, max_size = 7) {
  p + ggplot2::scale_size_area(name = "Retained samples", max_size = max_size, breaks = scales::breaks_pretty(n = 4)) + ggplot2::guides(size = ggplot2::guide_legend(order = 2, override.aes = list(alpha = 1)))
}

add_species_scales <- function(p, show_size = TRUE) {
  species_legend_rows <- dplyr::case_when(n_species_map <= 8 ~ 1L, n_species_map <= 18 ~ 2L, TRUE ~ 3L)
  species_guide <- ggplot2::guide_legend(order = 1, nrow = species_legend_rows, byrow = TRUE)
  p <- p + ggplot2::scale_colour_manual(name = "Species/group", values = species_label_palette, breaks = species_levels, labels = species_legend_labels, drop = FALSE, guide = species_guide) + ggplot2::scale_fill_manual(name = "Species/group", values = species_palette, breaks = species_levels, labels = species_legend_labels, drop = FALSE, guide = species_guide) + ggplot2::scale_shape_manual(name = "Species/group", values = species_shapes, breaks = species_levels, labels = species_legend_labels,
                                                                                                                                                                                                                                                                                                                                                                                                         drop = FALSE, guide = species_guide) + ggplot2::theme(legend.position = "bottom", legend.box = "vertical", legend.justification = "left")
  if (show_size) {
    p <- add_sample_size_scale(p, max_size = 7)
  }
  p
}

save_map_plot <- function(stem, plot, dimensions) {
  dimensions <- constrain_map_output_dimensions(dimensions)
  png_file <- file.path(maps_dir, paste0(stem, ".png"))
  ggplot2::ggsave(filename = png_file, plot = plot, width = unname(dimensions["width"]), height = unname(dimensions["height"]), units = "cm", dpi = map_dpi, bg = "white", limitsize = FALSE)
  saved_map_files <<- c(saved_map_files, png_file)
  if (isTRUE(map_save_pdf)) {
    pdf_file <- file.path(maps_dir, paste0(stem, ".pdf"))
    pdf_saved <- tryCatch({
      ggplot2::ggsave(filename = pdf_file, plot = plot, width = unname(dimensions["width"]), height = unname(dimensions["height"]), units = "cm", bg = "white", limitsize = FALSE)
      TRUE
    }, error = function(e) {
      warning("PNG was saved, but PDF failed for ", stem, ": ", conditionMessage(e), call. = FALSE)
      FALSE
    })
    if (pdf_saved)
      saved_map_files <<- c(saved_map_files, pdf_file)
  }
  invisible(png_file)
}

choose_facet_grid <- function(n_facets, max_columns = 5L) {
  max_columns <- max(1L, min(as.integer(max_columns), as.integer(n_facets)))
  candidate_columns <- seq_len(max_columns)
  candidate_rows <- ceiling(n_facets/candidate_columns)
  empty_panels <- candidate_columns * candidate_rows - n_facets
  target_grid_aspect <- 1.5
  grid_aspect <- candidate_columns/candidate_rows
  layout_score <- abs(log(grid_aspect/target_grid_aspect)) + 0.8 * empty_panels/n_facets + ifelse(candidate_columns < candidate_rows, 0.12, 0)
  best <- order(layout_score, empty_panels, abs(candidate_columns - candidate_rows), -candidate_columns)[1]
  c(ncol = candidate_columns[best], nrow = candidate_rows[best], empty = empty_panels[best])
}

map_inset_distance_matrix_km <- function(long, lat) {
  long <- suppressWarnings(as.numeric(long))
  lat <- suppressWarnings(as.numeric(lat))
  n <- length(long)
  output <- matrix(NA_real_, nrow = n, ncol = n)

  if (n == 0) return(output)

  earth_radius_km <- 6371.0088
  longitude_radians <- long * pi / 180
  latitude_radians <- lat * pi / 180

  for (index in seq_len(n)) {
    longitude_difference <- longitude_radians - longitude_radians[index]
    latitude_difference <- latitude_radians - latitude_radians[index]
    haversine_value <- sin(latitude_difference / 2)^2 +
      cos(latitude_radians[index]) * cos(latitude_radians) *
      sin(longitude_difference / 2)^2
    haversine_value <- pmin(1, pmax(0, haversine_value))

    output[index, ] <- 2 * earth_radius_km * atan2(
      sqrt(haversine_value),
      sqrt(1 - haversine_value)
    )
  }

  output
}

find_map_inset_clusters <- function(site_df) {
  coordinates <- site_df %>%
    dplyr::transmute(
      site_label = as.character(.data$site_label),
      long = suppressWarnings(as.numeric(.data$long)),
      lat = suppressWarnings(as.numeric(.data$lat))
    ) %>%
    dplyr::filter(is.finite(.data$long), is.finite(.data$lat)) %>%
    dplyr::distinct(.data$site_label, .keep_all = TRUE)

  minimum_sites <- max(2L, as.integer(map_inset_min_sites))

  if (nrow(coordinates) < minimum_sites) return(list())

  distance_km <- map_inset_distance_matrix_km(
    coordinates$long,
    coordinates$lat
  )
  finite_distances <- distance_km[is.finite(distance_km)]

  if (length(finite_distances) == 0) return(list())

  map_diameter_km <- max(finite_distances, na.rm = TRUE)
  radius_km <- suppressWarnings(as.numeric(map_inset_radius_km))[1]

  if (!is.finite(radius_km) || radius_km <= 0) {
    radius_km <- min(60, max(5, map_diameter_km * 0.05))
  }

  map_x_span <- max(diff(range(coordinates$long)), 0.20)
  map_y_span <- max(diff(range(coordinates$lat)), 0.20)
  x_relative <- (coordinates$long - min(coordinates$long)) / map_x_span
  y_relative <- (coordinates$lat - min(coordinates$lat)) / map_y_span
  visual_distance <- sqrt(
    outer(x_relative, x_relative, "-")^2 +
      outer(y_relative, y_relative, "-")^2
  )
  visual_trigger <- map_inset_numeric_setting(
    "map_inset_visual_trigger_prop",
    0.045
  )
  minimum_close_pairs <- max(
    1L,
    as.integer(map_inset_numeric_setting(
      "map_inset_min_close_pairs",
      3L
    ))
  )
  maximum_clusters <- max(
    0L,
    min(3L, as.integer(map_inset_numeric_setting(
      "map_inset_max_clusters",
      2L
    )))
  )
  separation_factor <- max(
    1,
    map_inset_numeric_setting(
      "map_inset_cluster_separation_factor",
      1.75
    )
  )
  halo_steps <- suppressWarnings(as.integer(map_inset_halo_steps))[1]
  if (!is.finite(halo_steps)) halo_steps <- 0L
  halo_steps <- max(0L, halo_steps)

  if (maximum_clusters == 0L) return(list())

  diagonal_distance <- distance_km
  diag(diagonal_distance) <- Inf
  available <- seq_len(nrow(coordinates))
  accepted <- list()
  accepted_seeds <- integer(0)
  attempts <- 0L

  while (length(available) >= minimum_sites &&
         length(accepted) < maximum_clusters &&
         attempts < nrow(coordinates) * 2L) {
    attempts <- attempts + 1L
    available_distance <- diagonal_distance[available, available, drop = FALSE]
    neighbour_count <- rowSums(
      available_distance <= radius_km,
      na.rm = TRUE
    ) + 1L
    seed_index <- available[which.max(neighbour_count)]
    cluster_index <- available[
      distance_km[seed_index, available] <= radius_km |
        available == seed_index
    ]

    if (halo_steps > 0L) {
      for (halo_step in seq_len(halo_steps)) {
        adjoining_index <- available[
          rowSums(
            distance_km[available, cluster_index, drop = FALSE] <= radius_km,
            na.rm = TRUE
          ) > 0
        ]
        expanded_index <- sort(unique(c(cluster_index, adjoining_index)))
        if (identical(expanded_index, cluster_index)) break
        cluster_index <- expanded_index
      }
    }

    if (length(cluster_index) < minimum_sites) {
      available <- setdiff(available, seed_index)
      next
    }

    cluster_visual_distance <- visual_distance[
      cluster_index,
      cluster_index,
      drop = FALSE
    ]
    close_links <- cluster_visual_distance <= visual_trigger
    diag(close_links) <- FALSE
    close_pair_count <- sum(close_links[upper.tri(close_links)], na.rm = TRUE)
    crowded_site_count <- sum(rowSums(close_links, na.rm = TRUE) > 0L)
    inset_is_needed <- crowded_site_count >= minimum_sites &&
      close_pair_count >= minimum_close_pairs

    sufficiently_separate <- TRUE
    if (length(accepted_seeds) > 0L) {
      sufficiently_separate <- all(
        distance_km[seed_index, accepted_seeds] >
          radius_km * separation_factor
      )
    }

    if (!inset_is_needed || !sufficiently_separate) {
      available <- if (sufficiently_separate) {
        setdiff(available, seed_index)
      } else {
        setdiff(available, cluster_index)
      }
      next
    }

    cluster_sites <- coordinates[cluster_index, , drop = FALSE]
    x_range <- range(cluster_sites$long)
    y_range <- range(cluster_sites$lat)
    x_span <- diff(x_range)
    y_span <- diff(y_range)
    mean_latitude <- mean(cluster_sites$lat)
    minimum_x_padding <- max(
      0.015,
      (radius_km /
         (111.32 * max(0.20, cos(mean_latitude * pi / 180)))) * 0.10
    )
    minimum_y_padding <- max(0.015, (radius_km / 110.57) * 0.10)
    x_padding <- max(x_span * map_inset_padding, minimum_x_padding)
    y_padding <- max(y_span * map_inset_padding, minimum_y_padding)

    accepted[[length(accepted) + 1L]] <- list(
      site_labels = cluster_sites$site_label,
      n_sites = nrow(cluster_sites),
      radius_km = radius_km,
      close_pairs = close_pair_count,
      crowded_sites = crowded_site_count,
      density_score = close_pair_count + crowded_site_count,
      centre_site = coordinates$site_label[seed_index],
      xlim = c(x_range[1] - x_padding, x_range[2] + x_padding),
      ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)
    )
    accepted_seeds <- c(accepted_seeds, seed_index)
    available <- setdiff(available, cluster_index)
  }

  if (length(accepted) > 1L) {
    density_scores <- vapply(accepted, `[[`, numeric(1), "density_score")
    accepted <- accepted[order(density_scores, decreasing = TRUE)]
    secondary_minimum_sites <- max(
      minimum_sites,
      as.integer(map_inset_numeric_setting(
        "map_inset_secondary_min_sites",
        6L
      ))
    )
    secondary_minimum_pairs <- max(
      minimum_close_pairs,
      as.integer(map_inset_numeric_setting(
        "map_inset_secondary_min_close_pairs",
        8L
      ))
    )
    secondary_density_ratio <- min(
      1,
      max(0, map_inset_numeric_setting(
        "map_inset_secondary_density_ratio",
        0.55
      ))
    )
    primary_density <- accepted[[1]]$density_score
    retained_secondary <- Filter(function(cluster) {
      cluster$crowded_sites >= secondary_minimum_sites &&
        cluster$close_pairs >= secondary_minimum_pairs &&
        cluster$density_score >= primary_density * secondary_density_ratio
    }, accepted[-1])
    accepted <- c(accepted[1], retained_secondary)
  }

  accepted
}

find_map_inset_cluster <- function(site_df) {
  clusters <- find_map_inset_clusters(site_df)
  if (length(clusters) == 0L) NULL else clusters[[1]]
}

map_inset_numeric_setting <- function(setting_name, default_value) {
  value <- get0(
    setting_name,
    envir = .GlobalEnv,
    inherits = TRUE,
    ifnotfound = default_value
  )
  value <- suppressWarnings(as.numeric(value))[1]
  if (!is.finite(value)) default_value else value
}

estimate_map_inset_land_fraction <- function(xmin,
                                             xmax,
                                             ymin,
                                             ymax,
                                             land_sf = NULL) {
  if (!requireNamespace("sf", quietly = TRUE) ||
      !all(is.finite(c(xmin, xmax, ymin, ymax))) ||
      xmax <= xmin || ymax <= ymin) {
    return(0.5)
  }

  if (is.null(land_sf)) {
    land_sf <- get0(
      "aus_states",
      envir = .GlobalEnv,
      inherits = TRUE,
      ifnotfound = NULL
    )
  }

  if (is.null(land_sf) || !inherits(land_sf, "sf")) return(0.5)

  grid_size <- max(
    6L,
    as.integer(map_inset_numeric_setting(
      "map_inset_land_grid_size",
      14L
    ))
  )
  x_step <- (xmax - xmin) / grid_size
  y_step <- (ymax - ymin) / grid_size
  grid_coordinates <- expand.grid(
    long = seq(
      xmin + x_step / 2,
      xmax - x_step / 2,
      length.out = grid_size
    ),
    lat = seq(
      ymin + y_step / 2,
      ymax - y_step / 2,
      length.out = grid_size
    )
  )

  tryCatch({
    land_sf <- suppressWarnings(sf::st_transform(land_sf, 4326))
    grid_sf <- sf::st_as_sf(
      grid_coordinates,
      coords = c("long", "lat"),
      crs = 4326,
      remove = FALSE
    )
    on_land <- lengths(suppressWarnings(
      sf::st_intersects(grid_sf, land_sf)
    )) > 0L
    mean(on_land)
  }, error = function(error_condition) {
    warning(
      "Inset ocean scoring could not evaluate the basemap; using neutral ",
      "land coverage for this run. Reason: ",
      conditionMessage(error_condition),
      call. = FALSE
    )
    0.5
  })
}

estimate_map_inset_land_fractions <- function(candidates,
                                              land_sf = NULL) {
  if (nrow(candidates) == 0L) return(numeric(0))
  if (!requireNamespace("sf", quietly = TRUE)) {
    return(rep(0.5, nrow(candidates)))
  }
  if (is.null(land_sf)) {
    land_sf <- get0(
      "aus_states",
      envir = .GlobalEnv,
      inherits = TRUE,
      ifnotfound = NULL
    )
  }
  if (is.null(land_sf) || !inherits(land_sf, "sf")) {
    return(rep(0.5, nrow(candidates)))
  }

  grid_size <- max(
    6L,
    as.integer(map_inset_numeric_setting("map_inset_land_grid_size", 14L))
  )
  tryCatch({
    land_sf <- suppressWarnings(sf::st_transform(land_sf, 4326))
    grids <- lapply(seq_len(nrow(candidates)), function(index) {
      x_step <- (candidates$xmax_data[index] -
        candidates$xmin_data[index]) / grid_size
      y_step <- (candidates$ymax_data[index] -
        candidates$ymin_data[index]) / grid_size
      expand.grid(
        candidate_id = index,
        long = seq(
          candidates$xmin_data[index] + x_step / 2,
          candidates$xmax_data[index] - x_step / 2,
          length.out = grid_size
        ),
        lat = seq(
          candidates$ymin_data[index] + y_step / 2,
          candidates$ymax_data[index] - y_step / 2,
          length.out = grid_size
        )
      )
    })
    grid_coordinates <- dplyr::bind_rows(grids)
    grid_sf <- sf::st_as_sf(
      grid_coordinates,
      coords = c("long", "lat"),
      crs = 4326,
      remove = FALSE
    )
    grid_coordinates$on_land <- lengths(suppressWarnings(
      sf::st_intersects(grid_sf, land_sf)
    )) > 0L
    fractions <- tapply(
      grid_coordinates$on_land,
      grid_coordinates$candidate_id,
      mean
    )
    as.numeric(fractions[as.character(seq_len(nrow(candidates)))])
  }, error = function(error_condition) {
    warning(
      "Batch inset ocean scoring failed; using individual estimates. ",
      "Reason: ", conditionMessage(error_condition),
      call. = FALSE
    )
    vapply(seq_len(nrow(candidates)), function(index) {
      estimate_map_inset_land_fraction(
        candidates$xmin_data[index],
        candidates$xmax_data[index],
        candidates$ymin_data[index],
        candidates$ymax_data[index],
        land_sf = land_sf
      )
    }, numeric(1))
  })
}

map_inset_rectangles_intersect_land <- function(candidates,
                                                land_sf = NULL) {
  if (!requireNamespace("sf", quietly = TRUE) || nrow(candidates) == 0L) {
    return(rep(FALSE, nrow(candidates)))
  }
  if (is.null(land_sf)) {
    land_sf <- get0(
      "aus_states",
      envir = .GlobalEnv,
      inherits = TRUE,
      ifnotfound = NULL
    )
  }
  if (is.null(land_sf) || !inherits(land_sf, "sf")) {
    return(rep(FALSE, nrow(candidates)))
  }

  tryCatch({
    land_sf <- suppressWarnings(sf::st_transform(land_sf, 4326))
    rectangles <- lapply(seq_len(nrow(candidates)), function(index) {
      sf::st_polygon(list(matrix(c(
        candidates$xmin_data[index], candidates$ymin_data[index],
        candidates$xmax_data[index], candidates$ymin_data[index],
        candidates$xmax_data[index], candidates$ymax_data[index],
        candidates$xmin_data[index], candidates$ymax_data[index],
        candidates$xmin_data[index], candidates$ymin_data[index]
      ), ncol = 2, byrow = TRUE)))
    })
    rectangle_sf <- sf::st_sf(
      geometry = sf::st_sfc(rectangles, crs = 4326)
    )
    lengths(suppressWarnings(sf::st_intersects(rectangle_sf, land_sf))) > 0L
  }, error = function(error_condition) {
    warning(
      "Exact inset land-intersection scoring failed; using sampled land ",
      "coverage instead. Reason: ",
      conditionMessage(error_condition),
      call. = FALSE
    )
    candidates$land_fraction > 0
  })
}

find_map_ocean_expansion_side <- function(limits, land_sf = NULL) {
  x_span <- diff(limits$xlim)
  y_span <- diff(limits$ylim)
  if (!all(is.finite(c(x_span, y_span))) || x_span <= 0 || y_span <= 0) {
    return(NULL)
  }
  strip_x <- max(0.05, x_span * 0.18)
  strip_y <- max(0.05, y_span * 0.18)
  side_boxes <- list(
    right = c(
      limits$xlim[2], limits$xlim[2] + strip_x,
      limits$ylim[1], limits$ylim[2]
    ),
    left = c(
      limits$xlim[1] - strip_x, limits$xlim[1],
      limits$ylim[1], limits$ylim[2]
    ),
    bottom = c(
      limits$xlim[1], limits$xlim[2],
      limits$ylim[1] - strip_y, limits$ylim[1]
    ),
    top = c(
      limits$xlim[1], limits$xlim[2],
      limits$ylim[2], limits$ylim[2] + strip_y
    )
  )
  water_fraction <- vapply(side_boxes, function(box) {
    1 - estimate_map_inset_land_fraction(
      box[1], box[2], box[3], box[4], land_sf = land_sf
    )
  }, numeric(1))
  ocean_threshold <- map_inset_numeric_setting(
    "map_inset_ocean_detection_threshold",
    0.65
  )
  best_side <- names(which.max(water_fraction))[1]
  if (length(best_side) == 0L ||
      !is.finite(water_fraction[best_side]) ||
      water_fraction[best_side] < ocean_threshold) {
    return(NULL)
  }
  list(
    side = best_side,
    water_fraction = unname(water_fraction[best_side]),
    all_water_fractions = water_fraction
  )
}

expand_map_limits_on_side <- function(limits, side, amount) {
  expanded <- limits
  amount <- max(0, suppressWarnings(as.numeric(amount))[1])
  if (!is.finite(amount) || amount <= 0) return(expanded)
  if (identical(side, "right")) {
    expanded$xlim[2] <- expanded$xlim[2] + amount
  } else if (identical(side, "left")) {
    expanded$xlim[1] <- expanded$xlim[1] - amount
  } else if (identical(side, "top")) {
    expanded$ylim[2] <- expanded$ylim[2] + amount
  } else if (identical(side, "bottom")) {
    expanded$ylim[1] <- expanded$ylim[1] - amount
  }
  expanded
}

expand_map_limits_toward_open_space <- function(limits) {
  minimum_aspect <- map_inset_numeric_setting(
    "map_output_min_aspect",
    2 / 3
  )
  maximum_aspect <- map_inset_numeric_setting(
    "map_output_max_aspect",
    3 / 2
  )
  longitude_scale <- max(0.20, cos(mean(limits$ylim) * pi / 180))
  geographic_aspect <-
    diff(limits$xlim) * longitude_scale / diff(limits$ylim)

  if (geographic_aspect < minimum_aspect) {
    required_x_span <- minimum_aspect * diff(limits$ylim) / longitude_scale
    extra_x <- max(0, required_x_span - diff(limits$xlim))
    left_land <- estimate_map_inset_land_fraction(
      limits$xlim[1] - extra_x,
      limits$xlim[1],
      limits$ylim[1],
      limits$ylim[2]
    )
    right_land <- estimate_map_inset_land_fraction(
      limits$xlim[2],
      limits$xlim[2] + extra_x,
      limits$ylim[1],
      limits$ylim[2]
    )
    horizontal_side <- if (right_land <= left_land) "right" else "left"
    return(constrain_map_geographic_limits(
      limits,
      minimum_aspect = minimum_aspect,
      maximum_aspect = maximum_aspect,
      horizontal_side = horizontal_side
    ))
  }

  if (geographic_aspect > maximum_aspect) {
    required_y_span <- diff(limits$xlim) * longitude_scale / maximum_aspect
    extra_y <- max(0, required_y_span - diff(limits$ylim))
    bottom_land <- estimate_map_inset_land_fraction(
      limits$xlim[1],
      limits$xlim[2],
      limits$ylim[1] - extra_y,
      limits$ylim[1]
    )
    top_land <- estimate_map_inset_land_fraction(
      limits$xlim[1],
      limits$xlim[2],
      limits$ylim[2],
      limits$ylim[2] + extra_y
    )
    vertical_side <- if (top_land <= bottom_land) "top" else "bottom"
    return(constrain_map_geographic_limits(
      limits,
      minimum_aspect = minimum_aspect,
      maximum_aspect = maximum_aspect,
      vertical_side = vertical_side
    ))
  }

  limits
}

map_inset_point_to_segment_distance <- function(point_x,
                                                point_y,
                                                start_x,
                                                start_y,
                                                end_x,
                                                end_y) {
  segment_x <- end_x - start_x
  segment_y <- end_y - start_y
  segment_length_squared <- segment_x^2 + segment_y^2

  if (!is.finite(segment_length_squared) || segment_length_squared <= 0) {
    return(sqrt((point_x - start_x)^2 + (point_y - start_y)^2))
  }

  projection <- (
    (point_x - start_x) * segment_x +
      (point_y - start_y) * segment_y
  ) / segment_length_squared
  projection <- pmin(1, pmax(0, projection))

  nearest_x <- start_x + projection * segment_x
  nearest_y <- start_y + projection * segment_y
  sqrt((point_x - nearest_x)^2 + (point_y - nearest_y)^2)
}

score_map_inset_candidates <- function(candidates,
                                       coordinates,
                                       map_limits,
                                       focus_xlim = NULL,
                                       focus_ylim = NULL,
                                       point_buffer = 0.018,
                                       land_sf = NULL,
                                       reserved_positions = list()) {
  x_span <- diff(map_limits$xlim)
  y_span <- diff(map_limits$ylim)

  candidates$xmin_data <- map_limits$xlim[1] + candidates$left * x_span
  candidates$xmax_data <- map_limits$xlim[1] + candidates$right * x_span
  candidates$ymin_data <- map_limits$ylim[1] + candidates$bottom * y_span
  candidates$ymax_data <- map_limits$ylim[1] + candidates$top * y_span

  coordinates$x_relative <-
    (coordinates$long - map_limits$xlim[1]) / x_span
  coordinates$y_relative <-
    (coordinates$lat - map_limits$ylim[1]) / y_span

  candidates$covered_sites <- vapply(
    seq_len(nrow(candidates)),
    function(index) {
      sum(
        coordinates$x_relative >= candidates$left[index] - point_buffer &
          coordinates$x_relative <= candidates$right[index] + point_buffer &
          coordinates$y_relative >= candidates$bottom[index] - point_buffer &
          coordinates$y_relative <= candidates$top[index] + point_buffer,
        na.rm = TRUE
      )
    },
    integer(1)
  )

  candidates$land_fraction <- estimate_map_inset_land_fractions(
    candidates,
    land_sf = land_sf
  )
  candidates$water_fraction <- 1 - candidates$land_fraction
  candidates$distance_to_focus <- 0
  candidates$connector_crossings <- 0L
  candidates$inset_overlaps <- 0L
  candidates$inset_spacing_violations <- 0L

  if (length(reserved_positions) > 0L) {
    inset_spacing <- map_inset_numeric_setting("map_inset_min_spacing", 0.065)
    candidates$inset_overlaps <- vapply(
      seq_len(nrow(candidates)),
      function(index) {
        sum(vapply(reserved_positions, function(reserved) {
          if (is.null(reserved) ||
              !all(c("left", "right", "bottom", "top") %in% names(reserved))) {
            return(FALSE)
          }
          horizontal_overlap <- min(
            candidates$right[index],
            reserved$right
          ) - max(candidates$left[index], reserved$left)
          vertical_overlap <- min(
            candidates$top[index],
            reserved$top
          ) - max(candidates$bottom[index], reserved$bottom)
          horizontal_overlap > -0.012 && vertical_overlap > -0.012
        }, logical(1)))
      },
      integer(1)
    )
    candidates$inset_spacing_violations <- vapply(
      seq_len(nrow(candidates)),
      function(index) {
        sum(vapply(reserved_positions, function(reserved) {
          if (is.null(reserved) ||
              !all(c("left", "right", "bottom", "top") %in% names(reserved))) {
            return(FALSE)
          }
          horizontal_overlap <- min(
            candidates$right[index],
            reserved$right
          ) - max(candidates$left[index], reserved$left)
          vertical_overlap <- min(
            candidates$top[index],
            reserved$top
          ) - max(candidates$bottom[index], reserved$bottom)
          horizontal_overlap > -inset_spacing &&
            vertical_overlap > -inset_spacing
        }, logical(1)))
      },
      integer(1)
    )
  }

  focus_is_valid <- !is.null(focus_xlim) && !is.null(focus_ylim) &&
    all(is.finite(c(focus_xlim, focus_ylim)))

  if (focus_is_valid) {
    focus_left <- (focus_xlim[1] - map_limits$xlim[1]) / x_span
    focus_right <- (focus_xlim[2] - map_limits$xlim[1]) / x_span
    focus_bottom <- (focus_ylim[1] - map_limits$ylim[1]) / y_span
    focus_top <- (focus_ylim[2] - map_limits$ylim[1]) / y_span
    focus_centre_x <- mean(c(focus_left, focus_right))
    focus_centre_y <- mean(c(focus_bottom, focus_top))
    connector_buffer <- map_inset_numeric_setting(
      "map_inset_connector_buffer",
      0.022
    )

    for (index in seq_len(nrow(candidates))) {
      inset_centre_x <- mean(c(
        candidates$left[index],
        candidates$right[index]
      ))
      inset_centre_y <- mean(c(
        candidates$bottom[index],
        candidates$top[index]
      ))
      candidates$distance_to_focus[index] <- sqrt(
        (inset_centre_x - focus_centre_x)^2 +
          (inset_centre_y - focus_centre_y)^2
      )

      source_x <- if (inset_centre_x < focus_centre_x) {
        focus_left
      } else {
        focus_right
      }
      target_x <- if (focus_centre_x < inset_centre_x) {
        candidates$left[index]
      } else {
        candidates$right[index]
      }

      outside_source <- !(
        coordinates$x_relative >= focus_left - point_buffer &
          coordinates$x_relative <= focus_right + point_buffer &
          coordinates$y_relative >= focus_bottom - point_buffer &
          coordinates$y_relative <= focus_top + point_buffer
      )
      outside_inset <- !(
        coordinates$x_relative >= candidates$left[index] - point_buffer &
          coordinates$x_relative <= candidates$right[index] + point_buffer &
          coordinates$y_relative >= candidates$bottom[index] - point_buffer &
          coordinates$y_relative <= candidates$top[index] + point_buffer
      )
      eligible <- outside_source & outside_inset

      top_distance <- map_inset_point_to_segment_distance(
        coordinates$x_relative,
        coordinates$y_relative,
        source_x,
        focus_top,
        target_x,
        candidates$top[index]
      )
      bottom_distance <- map_inset_point_to_segment_distance(
        coordinates$x_relative,
        coordinates$y_relative,
        source_x,
        focus_bottom,
        target_x,
        candidates$bottom[index]
      )

      candidates$connector_crossings[index] <- sum(
        eligible &
          pmin(top_distance, bottom_distance) <= connector_buffer,
        na.rm = TRUE
      )
    }
  }

  candidates$layout_score <-
    map_inset_numeric_setting("map_inset_site_penalty", 100) *
      candidates$covered_sites +
    map_inset_numeric_setting("map_inset_connector_penalty", 12) *
      candidates$connector_crossings +
    map_inset_numeric_setting("map_inset_land_penalty", 10) *
      candidates$land_fraction +
    map_inset_numeric_setting("map_inset_overlap_penalty", 1000) *
      candidates$inset_overlaps +
    map_inset_numeric_setting("map_inset_overlap_penalty", 1000) *
      candidates$inset_spacing_violations +
    map_inset_numeric_setting("map_inset_distance_penalty", 0.35) *
      candidates$distance_to_focus

  candidates
}

choose_map_inset_position <- function(point_df,
                                      map_limits,
                                      focus_xlim = NULL,
                                      focus_ylim = NULL,
                                      preferred_corner = NULL,
                                      inset_width = map_inset_width,
                                      inset_height = map_inset_height,
                                      margin = map_inset_outer_margin,
                                      point_buffer = 0.018,
                                      reserved_positions = list()) {
  candidates <- data.frame(
    corner = c(
      "top_left", "bottom_left", "top_right", "bottom_right",
      "middle_left", "middle_right", "top_middle", "bottom_middle"
    ),
    left = c(
      margin,
      margin,
      1 - margin - inset_width,
      1 - margin - inset_width,
      margin,
      1 - margin - inset_width,
      (1 - inset_width) / 2,
      (1 - inset_width) / 2
    ),
    bottom = c(
      1 - margin - inset_height,
      margin,
      1 - margin - inset_height,
      margin,
      (1 - inset_height) / 2,
      (1 - inset_height) / 2,
      1 - margin - inset_height,
      margin
    ),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      right = .data$left + inset_width,
      top = .data$bottom + inset_height
    )

  coordinates <- point_df %>%
    dplyr::transmute(
      long = suppressWarnings(as.numeric(.data$long)),
      lat = suppressWarnings(as.numeric(.data$lat))
    ) %>%
    dplyr::filter(is.finite(.data$long), is.finite(.data$lat)) %>%
    dplyr::distinct()

  candidates <- score_map_inset_candidates(
    candidates,
    coordinates,
    map_limits,
    focus_xlim = focus_xlim,
    focus_ylim = focus_ylim,
    point_buffer = point_buffer,
    reserved_positions = reserved_positions
  )

  if (!is.null(preferred_corner) &&
      preferred_corner %in% candidates$corner) {
    chosen <- candidates[
      candidates$corner == preferred_corner,
      ,
      drop = FALSE
    ]
  } else {
    chosen <- candidates[
      order(
        candidates$layout_score,
        candidates$covered_sites,
        candidates$connector_crossings,
        candidates$inset_overlaps,
        candidates$inset_spacing_violations,
        candidates$land_fraction,
        candidates$distance_to_focus
      )[1],
      ,
      drop = FALSE
    ]
  }

  as.list(chosen[1, ])
}

plan_map_inset_layout <- function(point_df, map_limits, dense_cluster) {
  initial_position <- choose_map_inset_position(
    point_df,
    map_limits,
    focus_xlim = dense_cluster$xlim,
    focus_ylim = dense_cluster$ylim
  )

  preferred_corner <- initial_position$corner
  horizontal_side <- if (grepl("left$", preferred_corner)) {
    "left"
  } else if (grepl("right$", preferred_corner)) {
    "right"
  } else {
    "middle"
  }
  vertical_side <- if (grepl("^top", preferred_corner)) {
    "top"
  } else if (grepl("^bottom", preferred_corner)) {
    "bottom"
  } else {
    "middle"
  }
  expansion_side <- if (horizontal_side != "middle") {
    horizontal_side
  } else {
    vertical_side
  }

  expanded_limits <- map_limits
  original_x_span <- diff(map_limits$xlim)
  original_y_span <- diff(map_limits$ylim)
  expanded_position <- NULL

  for (iteration in 0:30) {
    expanded_limits <- map_limits

    if (expansion_side == "right") {
      expanded_limits$xlim[2] <-
        map_limits$xlim[2] + original_x_span * 0.055 * iteration
    } else if (expansion_side == "left") {
      expanded_limits$xlim[1] <-
        map_limits$xlim[1] - original_x_span * 0.055 * iteration
    } else if (expansion_side == "top") {
      expanded_limits$ylim[2] <-
        map_limits$ylim[2] + original_y_span * 0.055 * iteration
    } else if (expansion_side == "bottom") {
      expanded_limits$ylim[1] <-
        map_limits$ylim[1] - original_y_span * 0.055 * iteration
    }

    expanded_position <- choose_map_inset_position(
      point_df,
      expanded_limits,
      focus_xlim = dense_cluster$xlim,
      focus_ylim = dense_cluster$ylim,
      preferred_corner = preferred_corner
    )

    if (expanded_position$covered_sites == 0) break
  }

  if (!is.null(expanded_position)) {
    extra_horizontal <- original_x_span * max(
      0,
      suppressWarnings(as.numeric(map_inset_extra_side_prop))[1]
    )
    extra_vertical <- original_y_span * max(
      0,
      suppressWarnings(as.numeric(map_inset_extra_vertical_prop))[1]
    )

    if (horizontal_side != "middle" &&
        is.finite(extra_horizontal) && extra_horizontal > 0) {
      if (horizontal_side == "right") {
        expanded_limits$xlim[2] <-
          expanded_limits$xlim[2] + extra_horizontal
      } else {
        expanded_limits$xlim[1] <-
          expanded_limits$xlim[1] - extra_horizontal
      }
    }

    if (vertical_side != "middle" &&
        is.finite(extra_vertical) && extra_vertical > 0) {
      if (vertical_side == "top") {
        expanded_limits$ylim[2] <-
          expanded_limits$ylim[2] + extra_vertical
      } else {
        expanded_limits$ylim[1] <-
          expanded_limits$ylim[1] - extra_vertical
      }
    }

    expanded_position <- choose_map_inset_position(
      point_df,
      expanded_limits,
      focus_xlim = dense_cluster$xlim,
      focus_ylim = dense_cluster$ylim,
      preferred_corner = preferred_corner
    )
  }

  if (is.null(expanded_position) || expanded_position$covered_sites > 0) {
    warning(
      "The labelled-map inset could not be cleared of all mapped sites.",
      call. = FALSE
    )
  }

  list(
    map_limits = expanded_limits,
    position = expanded_position,
    inset_corner = preferred_corner,
    horizontal_side = horizontal_side,
    vertical_side = vertical_side
  )
}

map_inset_position_sides <- function(position_name) {
  c(
    if (grepl("left$", position_name)) "left" else
      if (grepl("right$", position_name)) "right" else character(0),
    if (grepl("^top", position_name)) "top" else
      if (grepl("^bottom", position_name)) "bottom" else character(0)
  )
}

map_inset_dimensions_for_cluster <- function(cluster, map_limits) {
  mean_cluster_latitude <- mean(cluster$ylim)
  cluster_aspect <-
    diff(cluster$xlim) *
    max(0.20, cos(mean_cluster_latitude * pi / 180)) /
    diff(cluster$ylim)
  minimum_aspect <- map_inset_numeric_setting(
    "map_inset_min_frame_aspect",
    2 / 3
  )
  maximum_aspect <- map_inset_numeric_setting(
    "map_inset_max_frame_aspect",
    3 / 2
  )
  target_aspect <- min(
    maximum_aspect,
    max(minimum_aspect, cluster_aspect)
  )

  mean_map_latitude <- mean(map_limits$ylim)
  map_aspect <-
    diff(map_limits$xlim) *
    max(0.20, cos(mean_map_latitude * pi / 180)) /
    diff(map_limits$ylim)
  map_aspect <- max(map_aspect, 0.05)
  fraction_ratio <- target_aspect / map_aspect
  base_area <- map_inset_numeric_setting("map_inset_width", 0.30) *
    map_inset_numeric_setting("map_inset_height", 0.34)
  inset_width <- sqrt(base_area * fraction_ratio)
  inset_height <- sqrt(base_area / fraction_ratio)

  minimum_fraction <- map_inset_numeric_setting(
    "map_inset_min_fraction",
    0.18
  )
  maximum_fraction <- map_inset_numeric_setting(
    "map_inset_max_fraction",
    0.46
  )
  if (max(inset_width, inset_height) > maximum_fraction) {
    resize <- maximum_fraction / max(inset_width, inset_height)
    inset_width <- inset_width * resize
    inset_height <- inset_height * resize
  }
  if (min(inset_width, inset_height) < minimum_fraction) {
    resize <- minimum_fraction / min(inset_width, inset_height)
    inset_width <- inset_width * resize
    inset_height <- inset_height * resize
  }

  c(
    width = inset_width,
    height = inset_height,
    target_aspect = target_aspect
  )
}

plan_map_inset_side_stack <- function(point_df,
                                      map_limits,
                                      dense_clusters,
                                      preferred_side = NULL,
                                      inset_scale = 1) {
  inset_count <- length(dense_clusters)
  if (inset_count < 1L || inset_count > 3L) return(NULL)

  dimensions <- lapply(dense_clusters, function(cluster) {
    dimension <- map_inset_dimensions_for_cluster(cluster, map_limits)
    requested_scale <- max(
      1,
      suppressWarnings(as.numeric(inset_scale))[1]
    )
    maximum_fraction <- map_inset_numeric_setting(
      "map_inset_max_fraction",
      0.46
    )
    effective_scale <- min(
      requested_scale,
      maximum_fraction / max(dimension[c("width", "height")])
    )
    dimension[c("width", "height")] <-
      dimension[c("width", "height")] * effective_scale
    dimension
  })
  margin <- map_inset_numeric_setting("map_inset_outer_margin", 0.012)
  minimum_spacing <- map_inset_numeric_setting(
    "map_inset_min_spacing",
    0.065
  )
  coordinates <- point_df %>%
    dplyr::transmute(
      long = suppressWarnings(as.numeric(.data$long)),
      lat = suppressWarnings(as.numeric(.data$lat))
    ) %>%
    dplyr::filter(is.finite(.data$long), is.finite(.data$lat)) %>%
    dplyr::distinct()

  candidate_sets <- lapply(seq_along(dense_clusters), function(cluster_index) {
    inset_width <- dimensions[[cluster_index]]["width"]
    inset_height <- dimensions[[cluster_index]]["height"]
    centre_x <- seq(
      margin + inset_width / 2,
      1 - margin - inset_width / 2,
      length.out = 11L
    )
    centre_y <- seq(
      margin + inset_height / 2,
      1 - margin - inset_height / 2,
      length.out = 13L
    )
    candidate_centres <- expand.grid(
      centre_x = centre_x,
      centre_y = centre_y
    )

    if (!is.null(preferred_side)) {
      candidate_centres <- switch(
        preferred_side,
        right = dplyr::filter(candidate_centres, .data$centre_x >= 0.50),
        left = dplyr::filter(candidate_centres, .data$centre_x <= 0.50),
        top = dplyr::filter(candidate_centres, .data$centre_y >= 0.50),
        bottom = dplyr::filter(candidate_centres, .data$centre_y <= 0.50),
        candidate_centres
      )
    }

    candidates <- candidate_centres %>%
      dplyr::transmute(
        corner = paste0("free_", dplyr::row_number()),
        left = .data$centre_x - inset_width / 2,
        bottom = .data$centre_y - inset_height / 2,
        right = .data$centre_x + inset_width / 2,
        top = .data$centre_y + inset_height / 2
      )
    scored <- score_map_inset_candidates(
      candidates,
      coordinates,
      map_limits,
      focus_xlim = dense_clusters[[cluster_index]]$xlim,
      focus_ylim = dense_clusters[[cluster_index]]$ylim
    )
    scored$centre_x <- (scored$left + scored$right) / 2
    scored$centre_y <- (scored$bottom + scored$top) / 2
    coordinate_x <- (coordinates$long - map_limits$xlim[1]) /
      diff(map_limits$xlim)
    coordinate_y <- (coordinates$lat - map_limits$ylim[1]) /
      diff(map_limits$ylim)
    scored$nearest_site_distance <- vapply(
      seq_len(nrow(scored)),
      function(candidate_index) {
        min(sqrt(
          (coordinate_x - scored$centre_x[candidate_index])^2 +
            (coordinate_y - scored$centre_y[candidate_index])^2
        ))
      },
      numeric(1)
    )
    scored$land_intersection <- map_inset_rectangles_intersect_land(scored)
    scored
  })

  choose_joint_layout <- function(require_ocean = TRUE) {
    usable_sets <- lapply(candidate_sets, function(candidates) {
      usable <- candidates[candidates$covered_sites == 0L, , drop = FALSE]
      if (require_ocean) {
        usable <- usable[!usable$land_intersection, , drop = FALSE]
      }
      if (nrow(usable) == 0L) return(usable)
      usable <- usable[order(
        usable$connector_crossings,
        usable$land_intersection,
        usable$land_fraction,
        usable$distance_to_focus
      ), , drop = FALSE]
      candidate_limit <- if (inset_count >= 3L) 18L else 36L
      utils::head(usable, candidate_limit)
    })
    if (any(vapply(usable_sets, nrow, integer(1)) == 0L)) return(NULL)
    available_space_fraction <- mean(vapply(
      seq_along(usable_sets),
      function(cluster_index) {
        nrow(usable_sets[[cluster_index]]) /
          max(1L, nrow(candidate_sets[[cluster_index]]))
      },
      numeric(1)
    ))

    if (inset_count == 1L) {
      usable <- usable_sets[[1]]
      available_centre_x <- mean(usable$centre_x)
      available_centre_y <- mean(usable$centre_y)
      usable$available_centre_distance <- sqrt(
        (usable$centre_x - available_centre_x)^2 +
          (usable$centre_y - available_centre_y)^2
      )
      selected_index <- if (require_ocean) {
        order(
          usable$connector_crossings,
          usable$available_centre_distance,
          usable$distance_to_focus
        )[1]
      } else {
        order(
          -usable$nearest_site_distance,
          usable$connector_crossings,
          usable$available_centre_distance
        )[1]
      }
      position <- usable[selected_index, , drop = FALSE]
      horizontal_side <- if (position$centre_x >= 0.50) "right" else "left"
      vertical_side <- if (position$centre_y >= 0.50) "top" else "bottom"
      position$corner <- paste0(
        vertical_side,
        "_free_1_",
        horizontal_side
      )
      position <- as.list(position[1, ])
      side <- if (abs(position$centre_x - 0.5) >=
          abs(position$centre_y - 0.5)) {
        horizontal_side
      } else {
        vertical_side
      }
      return(list(
        side = side,
        positions = list(position),
        dimensions = dimensions,
        covered_sites = 0L,
        inset_overlaps = 0L,
        connector_crossings = position$connector_crossings,
        water_fraction = position$water_fraction,
        minimum_water_fraction = position$water_fraction,
        uses_ocean = require_ocean,
        available_space_fraction = available_space_fraction,
        inset_scale = inset_scale,
        total_score = position$layout_score
      ))
    }

    combinations <- do.call(
      expand.grid,
      c(lapply(usable_sets, function(candidates) seq_len(nrow(candidates))),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE)
    )
    if (nrow(combinations) == 0L) return(NULL)

    evaluate_combination <- function(combination_index) {
      positions <- lapply(seq_len(inset_count), function(cluster_index) {
        usable_sets[[cluster_index]][
          combinations[combination_index, cluster_index],
          ,
          drop = FALSE
        ]
      })
      pair_indices <- utils::combn(seq_len(inset_count), 2L)
      pair_dx <- numeric(ncol(pair_indices))
      pair_dy <- numeric(ncol(pair_indices))
      pair_distance <- numeric(ncol(pair_indices))

      for (pair_index in seq_len(ncol(pair_indices))) {
        first <- positions[[pair_indices[1, pair_index]]]
        second <- positions[[pair_indices[2, pair_index]]]
        horizontal_overlap <- min(first$right, second$right) -
          max(first$left, second$left)
        vertical_overlap <- min(first$top, second$top) -
          max(first$bottom, second$bottom)
        if (horizontal_overlap > -0.012 && vertical_overlap > -0.012) {
          return(c(valid = 0, score = Inf))
        }
        horizontal_gap <- max(0, -horizontal_overlap)
        vertical_gap <- max(0, -vertical_overlap)
        if (horizontal_gap < minimum_spacing &&
            vertical_gap < minimum_spacing) {
          return(c(valid = 0, score = Inf))
        }
        pair_dx[pair_index] <- abs(first$centre_x - second$centre_x)
        pair_dy[pair_index] <- abs(first$centre_y - second$centre_y)
        pair_distance[pair_index] <- sqrt(
          pair_dx[pair_index]^2 + pair_dy[pair_index]^2
        )
      }

      centres_x <- vapply(positions, function(position) {
        position$centre_x
      }, numeric(1))
      centres_y <- vapply(positions, function(position) {
        position$centre_y
      }, numeric(1))
      available_centres_x <- unlist(lapply(usable_sets, `[[`, "centre_x"))
      available_centres_y <- unlist(lapply(usable_sets, `[[`, "centre_y"))
      centroid_penalty <- sqrt(
        (mean(centres_x) - mean(available_centres_x))^2 +
          (mean(centres_y) - mean(available_centres_y))^2
      )
      diagonal_balance <- mean(abs(pair_dx - pair_dy))
      distance_balance <- if (length(pair_distance) > 1L) {
        stats::sd(pair_distance)
      } else {
        0
      }
      base_score <- sum(vapply(positions, function(position) {
        position$layout_score
      }, numeric(1))) - 4 * mean(vapply(positions, function(position) {
        position$nearest_site_distance
      }, numeric(1)))
      score <- base_score +
        22 * diagonal_balance +
        8 * centroid_penalty +
        5 * distance_balance -
        1.5 * mean(pair_distance)
      c(valid = 1, score = score)
    }

    evaluations <- vapply(
      seq_len(nrow(combinations)),
      evaluate_combination,
      numeric(2)
    )
    valid <- evaluations["valid", ] == 1
    if (!any(valid)) return(NULL)
    best_combination <- which(valid)[which.min(evaluations["score", valid])]
    positions <- lapply(seq_len(inset_count), function(cluster_index) {
      position <- usable_sets[[cluster_index]][
        combinations[best_combination, cluster_index],
        ,
        drop = FALSE
      ]
      centre_x <- position$centre_x
      centre_y <- position$centre_y
      horizontal_side <- if (centre_x >= 0.50) "right" else "left"
      vertical_side <- if (centre_y >= 0.50) "top" else "bottom"
      position$corner <- paste0(
        vertical_side,
        "_free_",
        cluster_index,
        "_",
        horizontal_side
      )
      as.list(position[1, ])
    })

    mean_x <- mean(vapply(positions, `[[`, numeric(1), "centre_x"))
    mean_y <- mean(vapply(positions, `[[`, numeric(1), "centre_y"))
    side <- if (abs(mean_x - 0.5) >= abs(mean_y - 0.5)) {
      if (mean_x >= 0.5) "right" else "left"
    } else {
      if (mean_y >= 0.5) "top" else "bottom"
    }
    list(
      side = side,
      positions = positions,
      dimensions = dimensions,
      covered_sites = 0L,
      inset_overlaps = 0L,
      connector_crossings = sum(vapply(
        positions,
        `[[`,
        integer(1),
        "connector_crossings"
      )),
      water_fraction = mean(vapply(
        positions,
        `[[`,
        numeric(1),
        "water_fraction"
      )),
      minimum_water_fraction = min(vapply(
        positions,
        `[[`,
        numeric(1),
        "water_fraction"
      )),
      uses_ocean = require_ocean,
      available_space_fraction = available_space_fraction,
      inset_scale = inset_scale,
      total_score = evaluations["score", best_combination]
    )
  }

  ocean_plan <- choose_joint_layout(require_ocean = TRUE)
  if (!is.null(ocean_plan)) return(ocean_plan)
  choose_joint_layout(require_ocean = FALSE)
}

plan_map_inset_layouts <- function(point_df,
                                   map_limits,
                                   dense_clusters) {
  if (length(dense_clusters) == 0L) {
    return(list(map_limits = map_limits, layouts = list()))
  }

  # Place every warranted inset using the same two-dimensional free-space
  # search. A single inset is centred within the available ocean region (or the
  # largest point-free gap); multiple insets are jointly distributed across
  # that space. The independent corner placement below remains a fallback.
  if (length(dense_clusters) >= 1L) {
    placement_limits <- map_limits
    stacked_plan <- plan_map_inset_side_stack(
      point_df,
      placement_limits,
      dense_clusters
    )

    # If adjacent ocean exists but the present map extent is too narrow to hold
    # the full inset rectangle, enlarge that ocean-facing edge and retry before
    # accepting any point-free inland placement.
    ocean_expansion <- find_map_ocean_expansion_side(map_limits)
    if (!is.null(ocean_expansion) &&
        (is.null(stacked_plan) || !isTRUE(stacked_plan$uses_ocean))) {
      expansion_step <- max(
        0.02,
        map_inset_numeric_setting("map_inset_ocean_expansion_step", 0.10)
      )
      expansion_attempts <- max(
        1L,
        as.integer(map_inset_numeric_setting(
          "map_inset_ocean_expansion_attempts",
          7L
        ))
      )
      base_expansion_span <- if (ocean_expansion$side %in%
          c("left", "right")) {
        diff(map_limits$xlim)
      } else {
        diff(map_limits$ylim)
      }
      for (attempt in seq_len(expansion_attempts)) {
        trial_limits <- expand_map_limits_on_side(
          map_limits,
          ocean_expansion$side,
          base_expansion_span * expansion_step * attempt
        )
        trial_limits <- constrain_map_geographic_limits(
          trial_limits,
          minimum_aspect = map_inset_numeric_setting(
            "map_output_min_aspect",
            2 / 3
          ),
          maximum_aspect = map_inset_numeric_setting(
            "map_output_max_aspect",
            3 / 2
          ),
          horizontal_side = if (ocean_expansion$side %in%
              c("left", "right")) ocean_expansion$side else "balanced",
          vertical_side = if (ocean_expansion$side %in%
              c("top", "bottom")) ocean_expansion$side else "balanced"
        )
        ocean_trial <- plan_map_inset_side_stack(
          point_df,
          trial_limits,
          dense_clusters,
          preferred_side = ocean_expansion$side
        )
        if (!is.null(ocean_trial) && isTRUE(ocean_trial$uses_ocean)) {
          placement_limits <- trial_limits
          stacked_plan <- ocean_trial
          break
        }
      }
    }

    if (!is.null(stacked_plan)) {
      maximum_adaptive_scale <- max(
        1,
        map_inset_numeric_setting("map_inset_max_adaptive_scale", 1.40)
      )
      blank_area_reference <- max(
        0.01,
        map_inset_numeric_setting("map_inset_blank_area_reference", 0.25)
      )
      adaptive_scale <- 1 + (maximum_adaptive_scale - 1) * min(
        1,
        stacked_plan$available_space_fraction / blank_area_reference
      )
      if (adaptive_scale > 1.01) {
        enlarged_plan <- plan_map_inset_side_stack(
          point_df,
          placement_limits,
          dense_clusters,
          inset_scale = adaptive_scale
        )
        if (!is.null(enlarged_plan) &&
            (!isTRUE(stacked_plan$uses_ocean) ||
              isTRUE(enlarged_plan$uses_ocean))) {
          stacked_plan <- enlarged_plan
        } else {
          adaptive_scale <- 1
        }
      }
      stacked_limits <- placement_limits
      original_x_span <- diff(map_limits$xlim)
      original_y_span <- diff(map_limits$ylim)
      extra_horizontal <- original_x_span * max(
        0,
        suppressWarnings(as.numeric(map_inset_extra_side_prop))[1]
      )
      extra_vertical <- original_y_span * max(
        0,
        suppressWarnings(as.numeric(map_inset_extra_vertical_prop))[1]
      )
      if (stacked_plan$side == "left" && is.finite(extra_horizontal)) {
        stacked_limits$xlim[1] <- stacked_limits$xlim[1] - extra_horizontal
      } else if (stacked_plan$side == "right" && is.finite(extra_horizontal)) {
        stacked_limits$xlim[2] <- stacked_limits$xlim[2] + extra_horizontal
      } else if (stacked_plan$side == "bottom" && is.finite(extra_vertical)) {
        stacked_limits$ylim[1] <- stacked_limits$ylim[1] - extra_vertical
      } else if (stacked_plan$side == "top" && is.finite(extra_vertical)) {
        stacked_limits$ylim[2] <- stacked_limits$ylim[2] + extra_vertical
      }
      stacked_limits <- constrain_map_geographic_limits(
        stacked_limits,
        minimum_aspect = map_inset_numeric_setting(
          "map_output_min_aspect",
          2 / 3
        ),
        maximum_aspect = map_inset_numeric_setting(
          "map_output_max_aspect",
          3 / 2
        ),
        horizontal_side = if (stacked_plan$side %in% c("left", "right")) {
          stacked_plan$side
        } else {
          "balanced"
        },
        vertical_side = if (stacked_plan$side %in% c("top", "bottom")) {
          stacked_plan$side
        } else {
          "balanced"
        }
      )

      final_stack <- plan_map_inset_side_stack(
        point_df,
        stacked_limits,
        dense_clusters,
        preferred_side = stacked_plan$side,
        inset_scale = adaptive_scale
      )
      if (is.null(final_stack) ||
          (isTRUE(stacked_plan$uses_ocean) &&
            !isTRUE(final_stack$uses_ocean))) {
        adaptive_scale <- 1
        final_stack <- plan_map_inset_side_stack(
          point_df,
          stacked_limits,
          dense_clusters,
          preferred_side = stacked_plan$side,
          inset_scale = adaptive_scale
        )
      }
      if (!is.null(final_stack)) {
        layouts <- lapply(seq_along(dense_clusters), function(cluster_index) {
          position <- final_stack$positions[[cluster_index]]
          position_sides <- map_inset_position_sides(position$corner)
          list(
            map_limits = stacked_limits,
            position = position,
            inset_dimensions = final_stack$dimensions[[cluster_index]],
            inset_scale = final_stack$inset_scale,
            inset_corner = position$corner,
            horizontal_side = intersect(
              position_sides,
              c("left", "right")
            )[1],
            vertical_side = intersect(
              position_sides,
              c("top", "bottom")
            )[1]
          )
        })
        return(list(map_limits = stacked_limits, layouts = layouts))
      }
    }
  }

  preferred_positions <- character(length(dense_clusters))
  reserved_positions <- list()

  for (cluster_index in seq_along(dense_clusters)) {
    cluster <- dense_clusters[[cluster_index]]
    inset_dimensions <- map_inset_dimensions_for_cluster(
      cluster,
      map_limits
    )
    selected_position <- choose_map_inset_position(
      point_df,
      map_limits,
      focus_xlim = cluster$xlim,
      focus_ylim = cluster$ylim,
      inset_width = inset_dimensions["width"],
      inset_height = inset_dimensions["height"],
      reserved_positions = reserved_positions
    )
    preferred_positions[cluster_index] <- selected_position$corner
    reserved_positions[[cluster_index]] <- selected_position
  }

  original_x_span <- diff(map_limits$xlim)
  original_y_span <- diff(map_limits$ylim)
  expansion_counts <- stats::setNames(integer(4), c(
    "left", "right", "top", "bottom"
  ))
  expanded_limits <- map_limits
  final_positions <- vector("list", length(dense_clusters))

  for (iteration in 0:30) {
    reserved_positions <- list()
    final_positions <- vector("list", length(dense_clusters))

    for (cluster_index in seq_along(dense_clusters)) {
      cluster <- dense_clusters[[cluster_index]]
      inset_dimensions <- map_inset_dimensions_for_cluster(
        cluster,
        expanded_limits
      )
      final_positions[[cluster_index]] <- choose_map_inset_position(
        point_df,
        expanded_limits,
        focus_xlim = cluster$xlim,
        focus_ylim = cluster$ylim,
        preferred_corner = preferred_positions[cluster_index],
        inset_width = inset_dimensions["width"],
        inset_height = inset_dimensions["height"],
        reserved_positions = reserved_positions
      )
      reserved_positions[[cluster_index]] <- final_positions[[cluster_index]]
    }

    covered <- vapply(
      final_positions,
      `[[`,
      integer(1),
      "covered_sites"
    )
    overlaps <- vapply(
      final_positions,
      `[[`,
      integer(1),
      "inset_overlaps"
    )

    if (all(covered == 0L) && all(overlaps == 0L)) break

    problem_indices <- which(covered > 0L | overlaps > 0L)
    sides_to_expand <- unique(unlist(lapply(
      preferred_positions[problem_indices],
      map_inset_position_sides
    )))
    if (length(sides_to_expand) == 0L) break
    expansion_counts[sides_to_expand] <-
      expansion_counts[sides_to_expand] + 1L

    expanded_limits <- map_limits
    expanded_limits$xlim[1] <- map_limits$xlim[1] -
      original_x_span * 0.055 * expansion_counts["left"]
    expanded_limits$xlim[2] <- map_limits$xlim[2] +
      original_x_span * 0.055 * expansion_counts["right"]
    expanded_limits$ylim[1] <- map_limits$ylim[1] -
      original_y_span * 0.055 * expansion_counts["bottom"]
    expanded_limits$ylim[2] <- map_limits$ylim[2] +
      original_y_span * 0.055 * expansion_counts["top"]
  }

  occupied_sides <- unique(unlist(lapply(
    preferred_positions,
    map_inset_position_sides
  )))
  extra_horizontal <- original_x_span * max(
    0,
    suppressWarnings(as.numeric(map_inset_extra_side_prop))[1]
  )
  extra_vertical <- original_y_span * max(
    0,
    suppressWarnings(as.numeric(map_inset_extra_vertical_prop))[1]
  )

  if ("left" %in% occupied_sides && is.finite(extra_horizontal)) {
    expanded_limits$xlim[1] <- expanded_limits$xlim[1] - extra_horizontal
  }
  if ("right" %in% occupied_sides && is.finite(extra_horizontal)) {
    expanded_limits$xlim[2] <- expanded_limits$xlim[2] + extra_horizontal
  }
  if ("bottom" %in% occupied_sides && is.finite(extra_vertical)) {
    expanded_limits$ylim[1] <- expanded_limits$ylim[1] - extra_vertical
  }
  if ("top" %in% occupied_sides && is.finite(extra_vertical)) {
    expanded_limits$ylim[2] <- expanded_limits$ylim[2] + extra_vertical
  }

  horizontal_constraint_side <- intersect(
    occupied_sides,
    c("left", "right")
  )
  vertical_constraint_side <- intersect(
    occupied_sides,
    c("top", "bottom")
  )
  expanded_limits <- constrain_map_geographic_limits(
    expanded_limits,
    minimum_aspect = map_inset_numeric_setting(
      "map_output_min_aspect",
      2 / 3
    ),
    maximum_aspect = map_inset_numeric_setting(
      "map_output_max_aspect",
      3 / 2
    ),
    horizontal_side = if (length(horizontal_constraint_side) == 1L) {
      horizontal_constraint_side
    } else {
      "balanced"
    },
    vertical_side = if (length(vertical_constraint_side) == 1L) {
      vertical_constraint_side
    } else {
      "balanced"
    }
  )

  reserved_positions <- list()
  layouts <- vector("list", length(dense_clusters))
  for (cluster_index in seq_along(dense_clusters)) {
    cluster <- dense_clusters[[cluster_index]]
    inset_dimensions <- map_inset_dimensions_for_cluster(
      cluster,
      expanded_limits
    )
    position <- choose_map_inset_position(
      point_df,
      expanded_limits,
      focus_xlim = cluster$xlim,
      focus_ylim = cluster$ylim,
      preferred_corner = preferred_positions[cluster_index],
      inset_width = inset_dimensions["width"],
      inset_height = inset_dimensions["height"],
      reserved_positions = reserved_positions
    )
    reserved_positions[[cluster_index]] <- position
    position_sides <- map_inset_position_sides(position$corner)
    layouts[[cluster_index]] <- list(
      map_limits = expanded_limits,
      position = position,
      inset_dimensions = inset_dimensions,
      inset_corner = position$corner,
      horizontal_side = intersect(position_sides, c("left", "right"))[1],
      vertical_side = intersect(position_sides, c("top", "bottom"))[1]
    )
  }

  unresolved <- vapply(layouts, function(layout) {
    layout$position$covered_sites > 0L ||
      layout$position$inset_overlaps > 0L
  }, logical(1))
  if (any(unresolved)) {
    warning(
      "One or more inset footprints could not be completely cleared of ",
      "mapped sites or the other inset.",
      call. = FALSE
    )
  }

  list(map_limits = expanded_limits, layouts = layouts)
}

fit_map_zoom_limits_to_inset <- function(focus_xlim,
                                         focus_ylim,
                                         inset_position) {
  focus_x_span <- diff(focus_xlim)
  focus_y_span <- diff(focus_ylim)
  frame_x_span <- inset_position$xmax_data - inset_position$xmin_data
  frame_y_span <- inset_position$ymax_data - inset_position$ymin_data

  all_spans <- c(
    focus_x_span,
    focus_y_span,
    frame_x_span,
    frame_y_span
  )

  if (!all(is.finite(all_spans)) || any(all_spans <= 0)) {
    return(list(xlim = focus_xlim, ylim = focus_ylim))
  }

  focus_latitude <- mean(focus_ylim)
  frame_latitude <- mean(c(
    inset_position$ymin_data,
    inset_position$ymax_data
  ))
  focus_longitude_scale <- max(
    0.20,
    cos(focus_latitude * pi / 180)
  )
  frame_longitude_scale <- max(
    0.20,
    cos(frame_latitude * pi / 180)
  )
  frame_aspect <-
    (frame_x_span * frame_longitude_scale) / frame_y_span
  focus_aspect <-
    (focus_x_span * focus_longitude_scale) / focus_y_span

  if (focus_aspect < frame_aspect) {
    required_x_span <-
      frame_aspect * focus_y_span / focus_longitude_scale
    x_midpoint <- mean(focus_xlim)
    focus_xlim <- x_midpoint + c(-0.5, 0.5) * required_x_span
  } else if (focus_aspect > frame_aspect) {
    required_y_span <-
      focus_x_span * focus_longitude_scale / frame_aspect
    y_midpoint <- mean(focus_ylim)
    focus_ylim <- y_midpoint + c(-0.5, 0.5) * required_y_span
  }

  list(xlim = focus_xlim, ylim = focus_ylim)
}

blend_map_zoom_limits <- function(original_limits,
                                  aspect_fitted_limits) {
  fit_strength <- suppressWarnings(as.numeric(
    map_inset_aspect_fit_strength
  ))[1]

  if (!is.finite(fit_strength)) fit_strength <- 1.00
  fit_strength <- min(1, max(0, fit_strength))

  blend_one_axis <- function(original_axis, fitted_axis) {
    midpoint <- mean(original_axis)
    original_span <- diff(original_axis)
    fitted_span <- diff(fitted_axis)
    blended_span <- original_span +
      fit_strength * (fitted_span - original_span)

    midpoint + c(-0.5, 0.5) * blended_span
  }

  list(
    xlim = blend_one_axis(
      original_limits$xlim,
      aspect_fitted_limits$xlim
    ),
    ylim = blend_one_axis(
      original_limits$ylim,
      aspect_fitted_limits$ylim
    )
  )
}

add_labelled_map_inset <- function(main_plot,
                                   zoom_plot,
                                   dense_cluster,
                                   inset_layout,
                                   zoom_limits) {
  inset_position <- inset_layout$position
  cluster_centre_x <- mean(zoom_limits$xlim)
  inset_centre_x <- mean(c(
    inset_position$xmin_data,
    inset_position$xmax_data
  ))
  source_x <- if (inset_centre_x < cluster_centre_x) {
    zoom_limits$xlim[1]
  } else {
    zoom_limits$xlim[2]
  }
  target_x <- if (cluster_centre_x < inset_centre_x) {
    inset_position$xmin_data
  } else {
    inset_position$xmax_data
  }

  source_box <- data.frame(
    xmin = zoom_limits$xlim[1],
    xmax = zoom_limits$xlim[2],
    ymin = zoom_limits$ylim[1],
    ymax = zoom_limits$ylim[2]
  )
  connector_lines <- data.frame(
    x = rep(source_x, 2),
    y = c(zoom_limits$ylim[2], zoom_limits$ylim[1]),
    xend = rep(target_x, 2),
    yend = c(
      inset_position$ymax_data,
      inset_position$ymin_data
    )
  )
  inset_frame <- data.frame(
    xmin = inset_position$xmin_data,
    xmax = inset_position$xmax_data,
    ymin = inset_position$ymin_data,
    ymax = inset_position$ymax_data
  )

  main_plot +
    ggplot2::geom_rect(
      data = source_box,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      inherit.aes = FALSE,
      fill = NA,
      colour = "#17202A",
      linewidth = 0.60,
      linetype = "22"
    ) +
    ggplot2::geom_segment(
      data = connector_lines,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        xend = .data$xend,
        yend = .data$yend
      ),
      inherit.aes = FALSE,
      colour = "#465661",
      linewidth = 0.42,
      lineend = "round"
    ) +
    ggplot2::annotation_custom(
      grob = ggplot2::ggplotGrob(zoom_plot),
      xmin = inset_position$xmin_data,
      xmax = inset_position$xmax_data,
      ymin = inset_position$ymin_data,
      ymax = inset_position$ymax_data
    ) +
    ggplot2::geom_rect(
      data = inset_frame,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      inherit.aes = FALSE,
      fill = NA,
      colour = "#263845",
      linewidth = 0.55
    )
}

make_labelled_map_dimensions <- function(limits,
                                         n_labels,
                                         minimum_width = 38) {
  dimensions <- map_dimensions(limits, n_labels = n_labels)

  if (unname(dimensions["width"]) < minimum_width) {
    scale_factor <- minimum_width / unname(dimensions["width"])
    dimensions <- c(
      width = minimum_width,
      height = min(
        82,
        unname(dimensions["height"]) * scale_factor
      )
    )
  }

  dimensions
}

save_labelled_map_plot <- function(stem, plot, dimensions) {
  dimensions <- constrain_map_output_dimensions(dimensions)
  png_file <- file.path(maps_dir, paste0(stem, ".png"))

  ggplot2::ggsave(
    filename = png_file,
    plot = plot,
    width = unname(dimensions["width"]),
    height = unname(dimensions["height"]),
    units = "cm",
    dpi = map_labelled_dpi,
    bg = "white",
    limitsize = FALSE
  )

  saved_map_files <<- c(saved_map_files, png_file)

  if (isTRUE(map_save_pdf)) {
    pdf_file <- file.path(maps_dir, paste0(stem, ".pdf"))

    pdf_saved <- tryCatch({
      ggplot2::ggsave(
        filename = pdf_file,
        plot = plot,
        width = unname(dimensions["width"]),
        height = unname(dimensions["height"]),
        units = "cm",
        bg = "white",
        limitsize = FALSE
      )
      TRUE
    }, error = function(error_condition) {
      warning(
        "PNG was saved, but PDF failed for ",
        stem,
        ": ",
        conditionMessage(error_condition),
        call. = FALSE
      )
      FALSE
    })

    if (pdf_saved) {
      saved_map_files <<- c(saved_map_files, pdf_file)
    }
  }

  invisible(png_file)
}

make_close_site_label_layer <- function(data,
                                        label_field,
                                        limits,
                                        text_size = 3,
                                        colour_field = NULL) {
  x_padding <- diff(limits$xlim) * 0.015
  y_padding <- diff(limits$ylim) * 0.015

  label_mapping <- if (is.null(colour_field)) {
    ggplot2::aes(
      x = .data$long,
      y = .data$lat,
      label = .data[[label_field]]
    )
  } else {
    ggplot2::aes(
      x = .data$long,
      y = .data$lat,
      label = .data[[label_field]],
      colour = .data[[colour_field]]
    )
  }

  label_args <- list(
    data = data,
    mapping = label_mapping,
    seed = map_seed,
    size = text_size,
    fontface = "bold",
    fill = "white",
    alpha = 1,
    label.padding = grid::unit(0.10, "lines"),
    box.padding = grid::unit(0.24, "lines"),
    point.padding = grid::unit(0.10, "lines"),
    force = 1.0,
    force_pull = 3.5,
    segment.linewidth = 0.24,
    min.segment.length = 0,
    max.time = 40,
    max.iter = 150000,
    max.overlaps = Inf,
    xlim = limits$xlim + c(1, -1) * x_padding,
    ylim = limits$ylim + c(1, -1) * y_padding,
    show.legend = FALSE
  )

  if (is.null(colour_field)) {
    label_args$colour <- "#263238"
    label_args$segment.colour <- "#607077"
  }

  do.call(ggrepel::geom_label_repel, label_args)
}

make_tight_labelled_map_limits <- function(site_df) {
  coordinates <- site_df %>%
    dplyr::transmute(
      long = suppressWarnings(as.numeric(.data$long)),
      lat = suppressWarnings(as.numeric(.data$lat))
    ) %>%
    dplyr::filter(is.finite(.data$long), is.finite(.data$lat))

  if (nrow(coordinates) == 0) {
    return(detail_limits)
  }

  x_range <- range(coordinates$long)
  y_range <- range(coordinates$lat)
  x_span <- max(diff(x_range), 0.20)
  y_span <- max(diff(y_range), 0.20)
  x_padding <- max(
    map_outer_minimum_padding,
    x_span * map_outer_padding_fraction
  )
  y_padding <- max(
    map_outer_minimum_padding,
    y_span * map_outer_padding_fraction
  )

  list(
    xlim = x_range + c(-x_padding, x_padding),
    ylim = y_range + c(-y_padding, y_padding)
  )
}

# Numbered-map helpers. Site numbers are stable within a run and always follow
# latitude from north to south (then longitude for ties), making dense maps
# readable without requiring every full site name to sit inside the map panel.
make_inset_number_label_layer <- function(data,
                                          zoom_limits,
                                          colour_field = NULL) {
  label_mapping <- if (is.null(colour_field)) {
    ggplot2::aes(
      x = .data$long,
      y = .data$lat,
      label = .data$site_number_label
    )
  } else {
    ggplot2::aes(
      x = .data$long,
      y = .data$lat,
      label = .data$site_number_label,
      colour = .data[[colour_field]]
    )
  }

  label_args <- list(
    data = data,
    mapping = label_mapping,
    seed = map_seed,
    size = map_number_label_text_size,
    fontface = "bold",
    fill = "white",
    alpha = 1,
    label.padding = grid::unit(0.10, "lines"),
    box.padding = grid::unit(0.18, "lines"),
    point.padding = grid::unit(0.08, "lines"),
    force = 1.0,
    force_pull = 3.5,
    segment.linewidth = 0.24,
    min.segment.length = 0,
    max.time = 25,
    max.iter = 100000,
    max.overlaps = Inf,
    xlim = zoom_limits$xlim,
    ylim = zoom_limits$ylim,
    show.legend = FALSE
  )
  if (is.null(colour_field)) {
    label_args$colour <- "#263238"
    label_args$segment.colour <- "#607077"
  }

  do.call(ggrepel::geom_label_repel, label_args)
}

make_numbered_site_key <- function(data, colour_field, palette, title) {
  key_data <- data %>%
    dplyr::transmute(
      site_number = as.integer(.data$site_number),
      key_y = -as.integer(.data$site_number),
      key_text = base::mapply(
        function(number, site, count) {
          name_lines <- base::strwrap(
            paste0(number, "  ", site),
            width = map_number_key_wrap_width
          )
          paste(
            c(
              name_lines,
              paste0("(n = ", format_sample_count(count), ")")
            ),
            collapse = "\n"
          )
        },
        .data$site_number,
        .data$display_site,
        .data$n_filtered,
        USE.NAMES = FALSE
      ),
      colour_group = as.character(.data[[colour_field]])
    ) %>%
    dplyr::arrange(.data$site_number)

  n_key_columns <- if (nrow(key_data) > 20L) 2L else 1L
  n_key_rows <- ceiling(nrow(key_data) / n_key_columns)
  key_column_spacing <- if (n_key_columns > 1L) {
    map_number_key_column_spacing
  } else {
    1
  }
  key_data <- key_data %>%
    dplyr::mutate(
      key_column = (dplyr::row_number() - 1L) %/% n_key_rows + 1L,
      key_row = (dplyr::row_number() - 1L) %% n_key_rows + 1L,
      key_x = (.data$key_column - 1L) * key_column_spacing,
      key_y = -.data$key_row
    )
  key_xmax <- if (n_key_columns > 1L) key_column_spacing + 1 else 1

  ggplot2::ggplot(key_data) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$key_x + 0.025,
        y = .data$key_y,
        fill = .data$colour_group,
        colour = .data$colour_group
      ),
      shape = 21,
      size = map_number_key_point_size,
      stroke = 0.40,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$key_x + 0.070,
        y = .data$key_y,
        label = .data$key_text
      ),
      hjust = 0,
      size = map_number_key_text_size,
      fontface = "bold",
      colour = "#263238",
      lineheight = 0.88
    ) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_colour_manual(values = palette, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(0, key_xmax), clip = "off") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0),
      plot.margin = ggplot2::margin(5.5, 3, 5.5, 0)
    )
}

make_numbered_site_map <- function(point_data,
                                   label_data,
                                   colour_field,
                                   palette,
                                   title,
                                   label_colour_field = colour_field) {
  main_label_data <- if (length(labelled_map_clusters) == 0L) {
    label_data
  } else {
    label_data %>%
      dplyr::filter(!.data$site_label %in% labelled_map_inset_sites)
  }

  output <- make_map_base(
    labelled_map_limits,
    title = title,
    subtitle = labelled_map_subtitle,
    caption = "Site numbers are ordered from northernmost to southernmost; see the key at right for full site names and retained sample counts.",
    projection_crs = 4326
  ) +
    ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(
        x = .data$long,
        y = .data$lat,
        fill = .data[[colour_field]],
        colour = .data[[colour_field]]
      ),
      shape = 21,
      size = map_labelled_point_size,
      stroke = 0.52,
      alpha = 1
    ) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_colour_manual(values = palette, guide = "none")

  if (length(labelled_map_clusters) > 0L) for (
    cluster_index in seq_along(labelled_map_clusters)
  ) {
    current_cluster <- labelled_map_clusters[[cluster_index]]
    current_layout <- labelled_map_layouts[[cluster_index]]
    current_zoom_limits <- labelled_map_zoom_limits_list[[cluster_index]]
    inset_points <- point_data %>%
      dplyr::filter(.data$site_label %in% current_cluster$site_labels)
    inset_labels <- label_data %>%
      dplyr::filter(.data$site_label %in% current_cluster$site_labels)

    zoom_map <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        data = aus_states,
        fill = "#F4F0E6",
        colour = "#6F7C82",
        linewidth = 0.28
      ) +
      ggplot2::geom_point(
        data = inset_points,
        ggplot2::aes(
          x = .data$long,
          y = .data$lat,
          fill = .data[[colour_field]],
          colour = .data[[colour_field]]
        ),
        shape = 21,
        size = map_labelled_point_size * map_inset_point_scale,
        stroke = 0.50,
        alpha = 1,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = palette, guide = "none") +
      ggplot2::scale_colour_manual(values = palette, guide = "none") +
      make_inset_number_label_layer(
        inset_labels,
        current_zoom_limits,
        colour_field = label_colour_field
      ) +
      ggplot2::coord_sf(
        crs = sf::st_crs(4326),
        default_crs = sf::st_crs(4326),
        datum = sf::st_crs(4326),
        xlim = current_zoom_limits$xlim,
        ylim = current_zoom_limits$ylim,
        expand = FALSE,
        clip = "on"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "#EAF3F5", colour = NA),
        plot.background = ggplot2::element_rect(fill = "#EAF3F5", colour = NA),
        plot.margin = ggplot2::margin(0, 0, 0, 0)
      )

    output <- add_labelled_map_inset(
      output,
      zoom_map,
      current_cluster,
      current_layout,
      current_zoom_limits
    )
  }

  output + make_close_site_label_layer(
    main_label_data,
    "site_number_label",
    limits = labelled_map_limits,
    text_size = map_number_label_text_size,
    colour_field = label_colour_field
  )
}


#### 3. Map code (evaluated by make_rnr_base_map() and run_rnr_maps()) ####
# Each chunk is the original MAPS-section code, unchanged. Chunks are run in
# rnr_maps_env, so the objects they create (map_df, site_centres,
# labelled_map_clusters, ...) are available to later workflow sections, as before.
rnr_map_chunks <- list()

# ---- Base map for later scatterpie maps (base_map, divxlims, divylims) ----
rnr_map_chunks$base_map <- quote({
# Required packages
if (!requireNamespace("sf", quietly = TRUE)) stop("The sf package is required for the map section. Install with install.packages('sf')", call. = FALSE)
if (!requireNamespace("ggnewscale", quietly = TRUE)) stop("The ggnewscale package is required for maps with independently coloured labels and points. Install with install.packages('ggnewscale')", call. = FALSE)

# Prepare mapping data
map_df <- filtered_site_summary %>%
  dplyr::filter(!is.na(long), !is.na(lat), is.finite(long), is.finite(lat)) %>%
  dplyr::mutate(long = as.numeric(long), lat = as.numeric(lat),
                site_label = as.character(.data[[site_col_name]]),
                species_label = as.character(.data[[species_col_name]]))
if (nrow(map_df) == 0) stop("No valid coordinates found in filtered_site_summary for mapping.", call. = FALSE)

# Label coordinates
site_label_df <- map_df %>%
  dplyr::group_by(site_label) %>%
  dplyr::summarise(long = mean(long, na.rm = TRUE), lat = mean(lat, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(lat))

species_label_df <- map_df %>%
  dplyr::group_by(site_label, species_label) %>%
  dplyr::summarise(long = mean(long, na.rm = TRUE), lat = mean(lat, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(lat))

# Site label colours
site_colour_lookup <- if (is.null(names(site_colours))) {
  stats::setNames(site_colours[seq_along(unique(map_df$site_label))], unique(map_df$site_label))
} else site_colours
site_label_df$site_colour_hex <- unname(site_colour_lookup[site_label_df$site_label])
site_label_df$site_colour_hex[is.na(site_label_df$site_colour_hex)] <- "black"
site_label_fill_lookup <- stats::setNames(scales::alpha(unname(site_colour_lookup), 0.15), names(site_colour_lookup))
site_label_df$site_fill_hex <- unname(site_label_fill_lookup[site_label_df$site_label])
site_label_df$site_fill_hex[is.na(site_label_df$site_fill_hex)] <- scales::alpha("white", 0.85)

# Species label colours
species_colour_lookup <- if (is.null(names(sp_colours))) {
  stats::setNames(sp_colours[seq_along(unique(map_df$species_label))], unique(map_df$species_label))
} else sp_colours
species_label_df$species_colour_hex <- unname(species_colour_lookup[species_label_df$species_label])
species_label_df$species_colour_hex[is.na(species_label_df$species_colour_hex)] <- "black"
species_label_fill_lookup <- stats::setNames(scales::alpha(unname(species_colour_lookup), 0.15), names(species_colour_lookup))
species_label_df$species_fill_hex <- unname(species_label_fill_lookup[species_label_df$species_label])
species_label_df$species_fill_hex[is.na(species_label_df$species_fill_hex)] <- scales::alpha("white", 0.85)

map_limits <- make_auto_map_limits(site_label_df, pad_prop = 0.055, min_pad_x = 0.08,
                                   min_pad_y = 0.08, target_aspect = 1.35)
divxlims <- map_limits$xlim; divylims <- map_limits$ylim

label_layout <- make_single_column_label_positions(site_label_df, divxlims, divylims, side = "right",
                                                   x_gap_prop = 0.1, top_bottom_pad_prop = 0.012,
                                                   blend_actual_lat = 0.12)
label_df <- label_layout$labels; divxlims <- label_layout$xlim; divylims <- label_layout$ylim
label_hjust <- label_layout$hjust

species_label_layout <- make_single_column_label_positions(species_label_df, map_limits$xlim, map_limits$ylim,
                                                           side = "right", x_gap_prop = 0.1,
                                                           top_bottom_pad_prop = 0.012,
                                                           blend_actual_lat = 0.12)
species_label_plot_df <- species_label_layout$labels
species_label_hjust <- species_label_layout$hjust
divxlims <- c(min(divxlims[1], species_label_layout$xlim[1]), max(divxlims[2], species_label_layout$xlim[2]))
divylims <- c(min(divylims[1], species_label_layout$ylim[1]), max(divylims[2], species_label_layout$ylim[2]))

# Automatic sizing
n_sites <- nrow(site_label_df); n_species <- length(unique(map_df$species_label))
max_label_chars <- max(nchar(as.character(site_label_df$site_label)), na.rm = TRUE)
pointsize <- dplyr::case_when(n_sites <= 30 ~ 3.2, n_sites <= 75 ~ 2.7, n_sites <= 150 ~ 2.1, TRUE ~ 1.5)
pointstroke <- dplyr::case_when(n_sites <= 30 ~ 0.85, n_sites <= 75 ~ 0.75, n_sites <= 150 ~ 0.65, TRUE ~ 0.55)
labelsize_base <- dplyr::case_when(n_sites <= 20 ~ 6.0, n_sites <= 35 ~ 5.6, n_sites <= 50 ~ 5.0,
                                   n_sites <= 75 ~ 4.4, n_sites <= 100 ~ 3.8, n_sites <= 150 ~ 3.2, TRUE ~ 2.6)
labelsize_name_adjustment <- dplyr::case_when(max_label_chars <= 12 ~ 1.00, max_label_chars <= 20 ~ 0.94,
                                              max_label_chars <= 30 ~ 0.86, max_label_chars <= 45 ~ 0.76, TRUE ~ 0.66)
label_density <- n_sites / diff(divylims)
labelsize_density_adjustment <- dplyr::case_when(label_density <= 10 ~ 1.00, label_density <= 16 ~ 0.94,
                                                 label_density <= 24 ~ 0.86, label_density <= 35 ~ 0.76, TRUE ~ 0.64)
labelsize <- max(1.7, min(labelsize_base * labelsize_name_adjustment * labelsize_density_adjustment, 5.2))
segment_lwd <- dplyr::case_when(n_sites <= 50 ~ 0.4, n_sites <= 100 ~ 0.35, TRUE ~ 0.20)
legend_cols_species <- dplyr::case_when(n_species <= 10 ~ 1, n_species <= 20 ~ 2, TRUE ~ 3)

# Base map
base_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = ozmaps::abs_ste, fill = "grey96", colour = "grey35", linewidth = 0.25) +
  ggplot2::coord_sf(xlim = divxlims, ylim = divylims, expand = FALSE) +
  ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.15),
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90, size = 6),
                 axis.text.y = ggplot2::element_text(size = 6), axis.title = ggplot2::element_blank(),
                 legend.background = ggplot2::element_blank(), legend.key = ggplot2::element_blank())
})

# ---- Shared preparation: package/input checks, coordinate QC, palettes, limits ----
rnr_map_chunks$prep <- quote({
map_packages <- c("ggplot2", "dplyr", "ggrepel", "sf", "ozmaps", "ggspatial", "colorspace", "scales")
missing_map_packages <- map_packages[!vapply(map_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_map_packages) > 0) {
  stop("Install the package(s) required by the map section: ", paste(missing_map_packages, collapse = ", "), call. = FALSE)
}
if (!exists("filtered_site_summary", inherits = TRUE)) {
  stop("`filtered_site_summary` must exist before the map section.", call. = FALSE)
}
if (!exists("output_dir", inherits = TRUE) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
  stop("`output_dir` must be a single valid folder path.", call. = FALSE)
}
map_input <- as.data.frame(filtered_site_summary, stringsAsFactors = FALSE)
required_map_columns <- c(site_col_name, species_col_name, "long", "lat")
missing_map_columns <- setdiff(required_map_columns, names(map_input))
if (length(missing_map_columns) > 0) {
  stop("`filtered_site_summary` is missing: ", paste(missing_map_columns, collapse = ", "), call. = FALSE)
}
maps_dir <- file.path(output_dir, "plots", "Maps")
dir.create(maps_dir, recursive = TRUE, showWarnings = FALSE)
old_map_manifest <- file.path(maps_dir, "Map_files_created.txt")
if (file.exists(old_map_manifest)) {
  file.remove(old_map_manifest)
}
map_input$.map_row <- seq_len(nrow(map_input))
map_input$site_label <- as.character(map_input[[site_col_name]])
map_input$species_label <- as.character(map_input[[species_col_name]])
map_input$site_label[is.na(map_input$site_label) | !nzchar(map_input$site_label)] <- "Unknown_site"
map_input$species_label[is.na(map_input$species_label) | !nzchar(map_input$species_label)] <- "Unknown_species"
map_input$long_raw <- map_input$long
map_input$lat_raw <- map_input$lat
map_input$long <- parse_map_coordinate(map_input$long_raw)
map_input$lat <- parse_map_coordinate(map_input$lat_raw)
map_input$coordinates_swapped <- is.finite(map_input$long) & is.finite(map_input$lat) & abs(map_input$long) >= 8 & abs(map_input$long) <= 56 & abs(map_input$lat) >= 95 & abs(map_input$lat) <= 170
if (any(map_input$coordinates_swapped)) {
  old_long <- map_input$long[map_input$coordinates_swapped]
  map_input$long[map_input$coordinates_swapped] <- map_input$lat[map_input$coordinates_swapped]
  map_input$lat[map_input$coordinates_swapped] <- old_long
}
map_input$latitude_sign_corrected <- is.finite(map_input$lat) & map_input$lat >= 8 & map_input$lat <= 56
map_input$lat[map_input$latitude_sign_corrected] <- -abs(map_input$lat[map_input$latitude_sign_corrected])
map_input$longitude_sign_corrected <- is.finite(map_input$long) & map_input$long <= -95 & map_input$long >= -170
map_input$long[map_input$longitude_sign_corrected] <- abs(map_input$long[map_input$longitude_sign_corrected])
map_input$valid_for_map <- is.finite(map_input$long) & is.finite(map_input$lat) & map_input$long >= 95 & map_input$long <= 170 & map_input$lat >= -56 & map_input$lat <= -8
map_input$coordinate_qc <- ifelse(map_input$valid_for_map, "included", "excluded: missing or outside broad Australian bounds")
map_input$coordinate_qc[map_input$valid_for_map & map_input$coordinates_swapped] <- "included: latitude/longitude swapped"
map_input$coordinate_qc[map_input$valid_for_map & map_input$latitude_sign_corrected & !map_input$coordinates_swapped] <- "included: latitude sign corrected"
map_input$coordinate_qc[map_input$valid_for_map & map_input$longitude_sign_corrected & !map_input$coordinates_swapped] <- "included: longitude sign corrected"
coordinate_qc_file <- file.path(maps_dir, "Map_coordinate_QC.csv")
utils::write.csv(map_input[, c(".map_row", species_col_name, site_col_name, "long_raw", "lat_raw", "long", "lat", "coordinates_swapped", "latitude_sign_corrected", "longitude_sign_corrected", "valid_for_map", "coordinate_qc")], coordinate_qc_file, row.names = FALSE, na = "")
n_corrected_coordinates <- sum(map_input$coordinates_swapped | map_input$latitude_sign_corrected | map_input$longitude_sign_corrected, na.rm = TRUE)
n_excluded_coordinates <- sum(!map_input$valid_for_map)
if (n_corrected_coordinates > 0) {
  message("Map coordinate QC corrected ", n_corrected_coordinates, " row(s); see ", coordinate_qc_file)
}
if (n_excluded_coordinates > 0) {
  warning("Map coordinate QC excluded ", n_excluded_coordinates, " row(s); see ", coordinate_qc_file, call. = FALSE)
}
map_df <- map_input[map_input$valid_for_map, , drop = FALSE]
if (nrow(map_df) == 0) {
  stop("No valid Australian coordinates remain after map coordinate QC. See: ", coordinate_qc_file, call. = FALSE)
}
if ("n_filtered" %in% names(map_df)) {
  map_df$n_filtered <- suppressWarnings(as.numeric(map_df$n_filtered))
} else {
  map_df$n_filtered <- 1
}
map_df$n_filtered[!is.finite(map_df$n_filtered) | map_df$n_filtered < 0] <- 1
map_df$plot_n <- pmax(map_df$n_filtered, 1)
map_df$display_site <- clean_map_label(map_df$site_label)
map_df$display_species <- clean_map_label(map_df$species_label)
site_centres <- dplyr::arrange(dplyr::summarise(dplyr::group_by(map_df, .data$site_label), display_site = dplyr::first(.data$display_site), long = stats::weighted.mean(.data$long, w = .data$plot_n, na.rm = TRUE), lat = stats::weighted.mean(.data$lat, w = .data$plot_n, na.rm = TRUE), n_filtered = sum(.data$n_filtered, na.rm = TRUE), plot_n = sum(.data$plot_n, na.rm = TRUE), n_species = dplyr::n_distinct(.data$species_label), .groups = "drop"), dplyr::desc(.data$lat), .data$long)
site_centres$display_site_n <- paste0(site_centres$display_site, " (", format_sample_count(site_centres$n_filtered), ")")
site_centres$display_site_count <- paste0(
  site_centres$display_site,
  " (",
  format_sample_count(site_centres$n_filtered),
  ")"
)
site_centres$display_site_n_local <- site_centres$display_site_n
species_site_labels <- dplyr::arrange(dplyr::summarise(dplyr::group_by(map_df, .data$site_label, .data$species_label), display_site = dplyr::first(.data$display_site), display_species = dplyr::first(.data$display_species), long = stats::weighted.mean(.data$long, w = .data$plot_n, na.rm = TRUE), lat = stats::weighted.mean(.data$lat, w = .data$plot_n, na.rm = TRUE), n_filtered = sum(.data$n_filtered, na.rm = TRUE), plot_n = sum(.data$plot_n, na.rm = TRUE), .groups = "drop"), dplyr::desc(.data$lat),
                                      .data$long, .data$species_label)
species_site_labels$display_species_site_n <- paste0(species_site_labels$display_site, " [", species_site_labels$display_species, "; n = ", format_sample_count(species_site_labels$n_filtered), "]")
species_site_labels$display_species_site_n_local <- species_site_labels$display_species_site_n
species_site_labels$display_site_n <- paste0(species_site_labels$display_site, " (", format_sample_count(species_site_labels$n_filtered), ")")
species_site_labels$display_site_count <- paste0(
  species_site_labels$display_site,
  " (",
  format_sample_count(species_site_labels$n_filtered),
  ")"
)
species_levels <- sort(unique(map_df$species_label))
site_levels <- site_centres$site_label
n_species_map <- length(species_levels)
n_sites_map <- nrow(site_centres)
n_populations_map <- nrow(species_site_labels)
site_palette <- complete_map_palette(site_levels, supplied = if (exists("site_colours", inherits = TRUE)) site_colours else NULL, chroma = 92, luminance = 54)
site_label_palette <- stats::setNames(colorspace::darken(unname(site_palette), amount = 0.12), names(site_palette))
species_palette <- complete_map_palette(species_levels, supplied = if (exists("sp_colours", inherits = TRUE)) sp_colours else NULL, chroma = 94, luminance = 52)
species_label_palette <- stats::setNames(colorspace::darken(unname(species_palette), amount = 0.12), names(species_palette))
species_shapes <- complete_map_shapes(species_levels, supplied = if (exists("sp_shapes", inherits = TRUE)) sp_shapes else NULL)
site_legend_labels <- stats::setNames(clean_map_label(site_levels), site_levels)
species_legend_labels <- stats::setNames(clean_map_label(species_levels), species_levels)
detail_limits <- make_map_limits(site_centres, pad_fraction = if (n_sites_map <= 30) 0.14 else 0.1, minimum_span = if (n_sites_map == 1) 0.5 else 0.36)
detail_limits <- ensure_minimum_map_aspect(detail_limits, minimum_aspect = map_min_detail_aspect)
full_australia_limits <- list(xlim = range(c(112, 154.5, detail_limits$xlim)), ylim = range(c(-44.3, -9.5, detail_limits$ylim)))
study_name <- if (exists("fullspecies", inherits = TRUE) && length(fullspecies) == 1L && !is.na(fullspecies) && nzchar(fullspecies)) {
  clean_map_label(fullspecies)
} else if (length(species_levels) == 1L) {
  clean_map_label(species_levels)
} else {
  "Genetic sampling"
}
map_data_note <- paste0(scales::comma(sum(map_df$n_filtered, na.rm = TRUE)), " retained samples; ", scales::comma(n_populations_map), " population", ifelse(n_populations_map == 1, "", "s"), "; ", scales::comma(n_sites_map), " site", ifelse(n_sites_map == 1, "", "s"), "; ", scales::comma(n_species_map), " species/group", ifelse(n_species_map == 1, "", "s"))
aus_states <- ozmaps::abs_ste
if (is.na(sf::st_crs(aus_states))) {
  sf::st_crs(aus_states) <- 4283
}
aus_states <- suppressWarnings(sf::st_transform(aus_states, 4326))

detail_dimensions <- map_dimensions(detail_limits, n_labels = n_sites_map)
if (unname(detail_dimensions["width"]) < map_min_detail_width_cm) {
  detail_scale <- map_min_detail_width_cm/unname(detail_dimensions["width"])
  detail_dimensions <- c(width = map_min_detail_width_cm, height = min(75, unname(detail_dimensions["height"]) * detail_scale))
}
context_dimensions <- c(width = 28, height = 23)
label_text_size <- dplyr::case_when(n_sites_map <= 20 ~ 3.8, n_sites_map <= 45 ~ 3.4, n_sites_map <= 80 ~ 3, n_sites_map <= 140 ~ 2.6, TRUE ~ 2.25)
saved_map_files <- character()
})

# ---- Local and ordered label layouts (Maps 01-04, 06, 07) ----
rnr_map_chunks$local_labels <- quote({
local_site_labels <- select_spatial_labels(site_centres, map_repel_label_cap)
local_species_labels <- select_spatial_labels(species_site_labels, map_repel_label_cap)
max_local_label_characters <- max(nchar(unlist(strsplit(c(local_site_labels$display_site_n_local, local_species_labels$display_species_site_n_local), "\n", fixed = TRUE))), na.rm = TRUE)
repel_x_inset <- min(0.16, max(0.05, 0.0025 * max_local_label_characters))
repel_y_inset <- 0.018
repel_label_xlim <- detail_limits$xlim + c(1, -1) * diff(detail_limits$xlim) * repel_x_inset
repel_label_ylim <- detail_limits$ylim + c(1, -1) * diff(detail_limits$ylim) * repel_y_inset
local_site_labels <- add_local_label_nudges(local_site_labels, detail_limits)
local_species_labels <- add_local_label_nudges(local_species_labels, detail_limits)
local_label_note <- if (nrow(local_site_labels) < n_sites_map) {
  paste0("Local labels show ", nrow(local_site_labels), " representative sites; ", "the ordered-label map lists all ", n_sites_map, " sites.")
} else {
  "Every site is labelled."
}
local_species_label_note <- if (nrow(local_species_labels) < n_populations_map) {
  paste0("Local labels show ", nrow(local_species_labels), " representative site × group populations; the ordered-label map lists all ", n_populations_map, ".")
} else {
  "Every site × group population is labelled separately."
}
species_label_text_size <- dplyr::case_when(nrow(local_species_labels) <= 35 ~ 3.4, nrow(local_species_labels) <= 70 ~ 2.8, nrow(local_species_labels) <= 110 ~ 2.45, TRUE ~ 2.15)
site_ordered_input <- site_centres
site_ordered_input$display_site <- site_ordered_input$display_site_n
ordered_layout <- make_ordered_label_column(site_ordered_input, detail_limits, label_text_size)
ordered_label_df <- ordered_layout$data
ordered_limits <- ordered_layout$limits
ordered_dimensions <- map_dimensions(ordered_limits, n_labels = n_sites_map, ordered_labels = TRUE)
species_ordered_input <- species_site_labels
species_ordered_input$display_site <- species_ordered_input$display_species_site_n
species_ordered_layout <- make_ordered_label_column(species_ordered_input, detail_limits, species_label_text_size)
species_ordered_label_df <- species_ordered_layout$data
species_ordered_limits <- species_ordered_layout$limits
species_ordered_dimensions <- map_dimensions(species_ordered_limits, n_labels = n_populations_map, ordered_labels = TRUE)
species_ordered_dimensions["height"] <- min(unname(species_ordered_dimensions["height"]), max(34, 14 + 0.3 * n_populations_map))
species_detail_dimensions <- map_dimensions(detail_limits, n_labels = nrow(local_species_labels))
minimum_species_width_cm <- max(44, map_min_detail_width_cm)
if (unname(species_detail_dimensions["width"]) < minimum_species_width_cm) {
  species_detail_scale <- minimum_species_width_cm/unname(species_detail_dimensions["width"])
  species_detail_dimensions <- c(width = minimum_species_width_cm, height = min(85, unname(species_detail_dimensions["height"]) * species_detail_scale))
}
})

# ---- Map 00: Australia context ----
rnr_map_chunks$map_00 <- quote({
context_map <- make_map_base(full_australia_limits, title = paste0(study_name, " — Australian sampling context"), subtitle = map_data_note, caption = paste0("The rectangle is the extent used by detailed maps. ", "Point area represents retained sample count.")) + ggplot2::geom_rect(data = data.frame(xmin = detail_limits$xlim[1], xmax = detail_limits$xlim[2], ymin = detail_limits$ylim[1], ymax = detail_limits$ylim[2]), ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
                                                                                                                                                                                                                                                                                          inherit.aes = FALSE, fill = scales::alpha("#D95F02", 0.06), colour = "#D95F02", linewidth = 0.65) + ggplot2::geom_point(data = site_centres, ggplot2::aes(x = .data$long, y = .data$lat, size = .data$plot_n), shape = 21, fill = scales::alpha("#D95F02", 0.78), colour = "white", stroke = 0.55, show.legend = TRUE)
context_map <- add_sample_size_scale(context_map, max_size = 6) + ggplot2::theme(legend.position = "bottom")
save_map_plot("Map_00_Australia_context", context_map, context_dimensions)
})

# ---- Map 01: sites with local labels ----
rnr_map_chunks$map_01 <- quote({
site_local_map <- make_map_base(detail_limits, title = paste0(study_name, " — sites with local labels"), subtitle = map_data_note, caption = paste0(local_label_note, " Label lines identify the true site coordinate; point area represents retained sample count.")) + ggplot2::geom_point(data = site_centres, ggplot2::aes(x = .data$long, y = .data$lat, colour = .data$site_label, fill = .data$site_label, size = .data$plot_n), shape = 21, alpha = 0.92, stroke = 0.65, show.legend = TRUE) + make_local_label_layer(local_site_labels,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "display_site_n_local", "site_label", label_text_size) + ggplot2::scale_colour_manual(values = site_label_palette, breaks = site_levels, labels = site_legend_labels, guide = "none") + ggplot2::scale_fill_manual(values = site_palette, breaks = site_levels, labels = site_legend_labels, guide = "none")
site_local_map <- add_sample_size_scale(site_local_map, max_size = 7) + ggplot2::theme(legend.position = "bottom")
save_map_plot("Map_01_Sites_local_labels", site_local_map, detail_dimensions)
})

# ---- Map 02: sites ordered north to south ----
rnr_map_chunks$map_02 <- quote({
site_ordered_map <- make_map_base(ordered_limits, title = paste0(study_name, " — all sites, ordered north to south"), subtitle = paste0(map_data_note, ". Labels are aligned for reliable reading."), caption = paste0("Every connector begins at the true site coordinate. ", "Point area represents retained sample count."), projection_crs = 4326) + ggplot2::geom_segment(data = ordered_label_df, ggplot2::aes(x = .data$long, y = .data$lat, xend = .data$label_long, yend = .data$label_lat, colour = .data$site_label),
                                                                                                                                                                                                                                                                                                                                                                               linewidth = 0.34, alpha = 0.78, lineend = "round") + ggplot2::geom_point(data = site_centres, ggplot2::aes(x = .data$long, y = .data$lat, colour = .data$site_label, fill = .data$site_label, size = .data$plot_n), shape = 21, alpha = 0.94, stroke = 0.65) + ggplot2::geom_point(data = ordered_label_df, ggplot2::aes(x = .data$label_long, y = .data$label_lat, colour = .data$site_label, fill = .data$site_label), shape = 21, size = 1.9, stroke = 0.45, show.legend = FALSE) + ggplot2::geom_label(data = ordered_label_df,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ggplot2::aes(x = .data$label_text_long, y = .data$label_lat, label = .data$display_site, colour = .data$site_label), hjust = 0, size = label_text_size, fill = scales::alpha("white", 0.94), linewidth = 0.18, label.padding = grid::unit(0.11, "lines"), fontface = "bold", show.legend = FALSE) + ggplot2::scale_colour_manual(values = site_label_palette, breaks = site_levels, labels = site_legend_labels, guide = "none") + ggplot2::scale_fill_manual(values = site_palette, breaks = site_levels, labels = site_legend_labels,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        guide = "none")
site_ordered_map <- add_sample_size_scale(site_ordered_map, max_size = 7) + ggplot2::theme(legend.position = "bottom")
save_map_plot("Map_02_Sites_ordered_labels", site_ordered_map, ordered_dimensions)
})

# ---- Map 03: species with local labels ----
rnr_map_chunks$map_03 <- quote({
species_local_map <- make_map_base(detail_limits, title = paste0(study_name, " — species/groups with local site labels"), subtitle = map_data_note, caption = paste0(local_species_label_note, " Labels and connectors match species/group colour; point area represents retained sample count.")) + ggplot2::geom_point(data = map_df, ggplot2::aes(x = .data$long, y = .data$lat, colour = .data$species_label, fill = .data$species_label, shape = .data$species_label, size = .data$plot_n), alpha = 0.9, stroke = 0.7) +
  make_local_label_layer(local_species_labels, "display_species_site_n_local", "species_label", species_label_text_size)
species_local_map <- add_species_scales(species_local_map, show_size = TRUE)
save_map_plot("Map_03_Species_local_labels", species_local_map, species_detail_dimensions)
})

# ---- Map 04: species ordered labels ----
rnr_map_chunks$map_04 <- quote({
species_ordered_map <- make_map_base(species_ordered_limits, title = paste0(study_name, " — species/groups and all populations"), subtitle = paste0(map_data_note, ". Each site × group population has its own north-to-south label."), caption = paste0("Connectors lead from true site coordinates to the aligned labels. ", "Labels and connectors match species/group colour."), projection_crs = 4326) + ggplot2::geom_segment(data = species_ordered_label_df, ggplot2::aes(x = .data$long, y = .data$lat, xend = .data$label_long,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  yend = .data$label_lat, colour = .data$species_label), linewidth = 0.34, alpha = 0.74, lineend = "round") + ggplot2::geom_point(data = map_df, ggplot2::aes(x = .data$long, y = .data$lat, colour = .data$species_label, fill = .data$species_label, shape = .data$species_label, size = .data$plot_n), alpha = 0.92, stroke = 0.7) + ggplot2::geom_point(data = species_ordered_label_df, ggplot2::aes(x = .data$label_long, y = .data$label_lat, colour = .data$species_label, fill = .data$species_label), shape = 21,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            size = 1.7, stroke = 0.4, show.legend = FALSE) + ggplot2::geom_label(data = species_ordered_label_df, ggplot2::aes(x = .data$label_text_long, y = .data$label_lat, label = .data$display_species_site_n, colour = .data$species_label), hjust = 0, size = species_label_text_size, fill = scales::alpha("white", 0.94), linewidth = 0.18, label.padding = grid::unit(0.11, "lines"), fontface = "bold", show.legend = FALSE)
species_ordered_map <- add_species_scales(species_ordered_map, show_size = TRUE)
save_map_plot("Map_04_Species_ordered_labels", species_ordered_map, species_ordered_dimensions)
})

# ---- Map 05: species without labels ----
rnr_map_chunks$map_05 <- quote({
species_no_labels_map <- make_map_base(detail_limits, title = paste0(study_name, " — uncluttered species/group distribution"), subtitle = map_data_note, caption = paste0("Points use exact mean population coordinates (no positional jitter). ", "Transparency helps reveal dense or overlapping locations.")) + ggplot2::geom_point(data = map_df, ggplot2::aes(x = .data$long, y = .data$lat, colour = .data$species_label, fill = .data$species_label, shape = .data$species_label, size = .data$plot_n), alpha = 0.74,
                                                                                                                                                                                                                                                                                                                                       stroke = 0.65)
species_no_labels_map <- add_species_scales(species_no_labels_map, show_size = TRUE)
save_map_plot("Map_05_Species_no_labels", species_no_labels_map, detail_dimensions)
})

# ---- Map 06: sampling effort ----
rnr_map_chunks$map_06 <- quote({
sampling_effort_map <- make_map_base(detail_limits, title = paste0(study_name, " — sampling effort"), subtitle = map_data_note, caption = paste0(local_label_note, " Both point area and colour represent the retained sample count at each site.")) + ggplot2::geom_point(data = site_centres, ggplot2::aes(x = .data$long, y = .data$lat, size = .data$plot_n, fill = .data$plot_n), shape = 21, colour = "white", alpha = 0.88, stroke = 0.65) + make_local_label_layer(local_site_labels, "display_site_n_local",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                           text_size = label_text_size) + ggplot2::scale_fill_viridis_c(name = "Retained samples", option = "C", begin = 0.1, end = 0.92)
sampling_effort_map <- add_sample_size_scale(sampling_effort_map, max_size = 9) + ggplot2::theme(legend.position = "bottom")
save_map_plot("Map_06_Sampling_effort", sampling_effort_map, detail_dimensions)
})

# ---- Map 07: species facets ----
rnr_map_chunks$map_07 <- quote({
if (n_species_map > 1) {
  facet_label_list <- lapply(split(species_site_labels, species_site_labels$species_label), select_spatial_labels, cap = map_facet_label_cap)
  facet_label_df <- dplyr::bind_rows(facet_label_list)
  facet_labeller <- stats::setNames(clean_map_label(species_levels), species_levels)
  facet_layout <- choose_facet_grid(n_species_map, max_columns = map_facet_max_columns)
  facet_ncol <- unname(facet_layout["ncol"])
  facet_nrow <- unname(facet_layout["nrow"])
  species_facet_map <- make_map_base(detail_limits, title = paste0(study_name, " — distribution by species/group"), subtitle = paste0("All panels use the same extent. Up to ", map_facet_label_cap, " geographically distributed sites are labelled per panel."), caption = paste0("Common limits make longitudinal and latitudinal separation directly comparable. ", "Point area represents retained sample count."), ornaments = FALSE) + ggplot2::geom_point(data = map_df, ggplot2::aes(x = .data$long, y = .data$lat,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              colour = .data$species_label, fill = .data$species_label, size = .data$plot_n), shape = 21, alpha = 0.84, stroke = 0.6, show.legend = FALSE) + ggrepel::geom_text_repel(data = facet_label_df, ggplot2::aes(x = .data$long, y = .data$lat, label = .data$display_site_n, colour = .data$species_label), seed = map_seed, size = 2.25, fontface = "bold", box.padding = grid::unit(0.42, "lines"), point.padding = grid::unit(0.22, "lines"), min.segment.length = 0, segment.alpha = 0.8, segment.size = 0.32,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      force = 1.85, force_pull = 0.18, max.time = 15, max.iter = 60000, max.overlaps = Inf, xlim = repel_label_xlim, ylim = repel_label_ylim, show.legend = FALSE) + ggplot2::scale_colour_manual(values = species_label_palette, guide = "none") + ggplot2::scale_fill_manual(values = species_palette, guide = "none") + ggplot2::scale_size_area(max_size = 5.5, guide = "none") + ggplot2::facet_wrap(ggplot2::vars(species_label), labeller = ggplot2::as_labeller(facet_labeller), ncol = facet_ncol) + ggplot2::theme(strip.background = ggplot2::element_rect(fill = "#31454F",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      colour = NA), strip.text = ggplot2::element_text(colour = "white", face = "bold", size = 8), axis.text = ggplot2::element_text(size = 6), axis.title = ggplot2::element_text(size = 7), panel.spacing = grid::unit(0.35, "cm"))
  facet_dimensions <- c(width = min(62, max(36, 14 * facet_ncol + 6)), height = min(75, max(28, 18 * facet_nrow + 8)))
  save_map_plot("Map_07_Species_facets", species_facet_map, facet_dimensions)
}
})

# ---- Shared sizes for Maps 08-09 ----
rnr_map_chunks$profile <- quote({
profile_axis_text_size <- dplyr::case_when(n_sites_map <= 35 ~ 7.5, n_sites_map <= 70 ~ 6.2, n_sites_map <= 120 ~ 5, TRUE ~ 4)
profile_dimensions <- c(width = min(48, max(28, 22 + 0.28 * max(nchar(site_centres$display_site_n)))), height = min(90, max(18, 8 + 0.48 * n_sites_map)))
})

# ---- Map 08: longitude profile ----
rnr_map_chunks$map_08 <- quote({
site_order_north_to_south <- site_centres$site_label
profile_df <- map_df
profile_df$site_axis <- factor(profile_df$site_label, levels = rev(site_order_north_to_south))
site_axis_labels <- stats::setNames(paste0(site_centres$display_site, "  (n = ", format_sample_count(site_centres$n_filtered), "; ", format(round(abs(site_centres$lat), 2), trim = TRUE), "°S)"), site_centres$site_label)
site_longitude_ranges <- dplyr::summarise(dplyr::group_by(profile_df, .data$site_label, .data$site_axis), long_min = min(.data$long, na.rm = TRUE), long_max = max(.data$long, na.rm = TRUE), .groups = "drop")
longitude_profile <- ggplot2::ggplot() + ggplot2::geom_segment(data = site_longitude_ranges, ggplot2::aes(x = .data$long_min, xend = .data$long_max, y = .data$site_axis, yend = .data$site_axis), colour = "#BCC7CB", linewidth = 0.65, lineend = "round") + ggplot2::geom_point(data = profile_df, ggplot2::aes(x = .data$long, y = .data$site_axis, colour = .data$species_label, fill = .data$species_label, shape = .data$species_label, size = .data$plot_n, group = .data$species_label), position = ggplot2::position_dodge(width = 0.56),
                                                                                                                                                                                                                                                                                  alpha = 0.9, stroke = 0.7) + ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 7), labels = format_longitude, expand = ggplot2::expansion(mult = c(0.04, 0.08))) + ggplot2::scale_y_discrete(labels = site_axis_labels) + ggplot2::labs(title = paste0(study_name, " — longitude profile"), subtitle = paste0("Sites are ordered north-to-south; horizontal position is true longitude. ", "Vertical spacing is rank, not geographic distance."), x = "Longitude", y = "Site (retained samples; mean latitude)",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       caption = paste0("This view is designed to expose longitudinal divisions hidden by a single label column. ", "Site labels give total retained samples; point area represents the retained count for each species/group.")) + ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "#E4EAEC", linewidth = 0.4), panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_line(colour = "#D5DEE1", linewidth = 0.4), axis.text.y = ggplot2::element_text(size = profile_axis_text_size,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               colour = "#2E3D43"), axis.text.x = ggplot2::element_text(size = 8, colour = "#48565C"), axis.title = ggplot2::element_text(face = "bold", colour = "#334148"), plot.title = ggplot2::element_text(colour = "#17252B", size = 15, face = "bold"), plot.subtitle = ggplot2::element_text(colour = "#48565C", size = 9, margin = ggplot2::margin(b = 6)), plot.caption = ggplot2::element_text(colour = "#65737A", size = 7, hjust = 0), legend.position = "bottom", legend.box = "vertical", plot.background = ggplot2::element_rect(fill = "white",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  colour = NA), plot.margin = ggplot2::margin(8, 10, 8, 8))
longitude_profile <- add_species_scales(longitude_profile, show_size = TRUE)
save_map_plot("Map_08_Longitude_profile", longitude_profile, profile_dimensions)
})

# ---- Map 09: latitude profile ----
rnr_map_chunks$map_09 <- quote({
site_centres_west_to_east <- dplyr::arrange(site_centres, .data$long, dplyr::desc(.data$lat))
site_order_west_to_east <- site_centres_west_to_east$site_label
latitude_profile_df <- map_df
latitude_profile_df$site_axis <- factor(latitude_profile_df$site_label, levels = rev(site_order_west_to_east))
latitude_site_axis_labels <- stats::setNames(paste0(site_centres_west_to_east$display_site, "  (n = ", format_sample_count(site_centres_west_to_east$n_filtered), "; ", format_longitude(site_centres_west_to_east$long), ")"), site_centres_west_to_east$site_label)
site_latitude_ranges <- dplyr::summarise(dplyr::group_by(latitude_profile_df, .data$site_label, .data$site_axis), lat_min = min(.data$lat, na.rm = TRUE), lat_max = max(.data$lat, na.rm = TRUE), .groups = "drop")
latitude_profile <- ggplot2::ggplot() + ggplot2::geom_segment(data = site_latitude_ranges, ggplot2::aes(x = .data$lat_min, xend = .data$lat_max, y = .data$site_axis, yend = .data$site_axis), colour = "#BCC7CB", linewidth = 0.65, lineend = "round") + ggplot2::geom_point(data = latitude_profile_df, ggplot2::aes(x = .data$lat, y = .data$site_axis, colour = .data$species_label, fill = .data$species_label, shape = .data$species_label, size = .data$plot_n, group = .data$species_label), position = ggplot2::position_dodge(width = 0.56),
                                                                                                                                                                                                                                                                              alpha = 0.9, stroke = 0.7) + ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 7), labels = format_latitude, expand = ggplot2::expansion(mult = c(0.04, 0.08))) + ggplot2::scale_y_discrete(labels = latitude_site_axis_labels) + ggplot2::labs(title = paste0(study_name, " — latitude profile"), subtitle = paste0("Sites are ordered west-to-east; horizontal position is true latitude. ", "Vertical spacing is rank, not geographic distance."), x = "Latitude", y = "Site (retained samples; mean longitude)",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           caption = paste0("This view is designed to expose latitudinal divisions hidden by a single label column. ", "Site labels give total retained samples; point area represents the retained count for each species/group.")) + ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "#E4EAEC", linewidth = 0.4), panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_line(colour = "#D5DEE1", linewidth = 0.4), axis.text.y = ggplot2::element_text(size = profile_axis_text_size,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  colour = "#2E3D43"), axis.text.x = ggplot2::element_text(size = 8, colour = "#48565C"), axis.title = ggplot2::element_text(face = "bold", colour = "#334148"), plot.title = ggplot2::element_text(colour = "#17252B", size = 15, face = "bold"), plot.subtitle = ggplot2::element_text(colour = "#48565C", size = 9, margin = ggplot2::margin(b = 6)), plot.caption = ggplot2::element_text(colour = "#65737A", size = 7, hjust = 0), legend.position = "bottom", legend.box = "vertical", plot.background = ggplot2::element_rect(fill = "white",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     colour = NA), plot.margin = ggplot2::margin(8, 10, 8, 8))
latitude_profile <- add_species_scales(latitude_profile, show_size = TRUE)
save_map_plot("Map_09_Latitude_profile", latitude_profile, profile_dimensions)
})

# ---- Shared inset plan (Maps 10-14) ----
rnr_map_chunks$inset <- quote({
#  Shared inset plan -

labelled_map_base_limits <- make_tight_labelled_map_limits(site_centres)
labelled_map_base_limits <- expand_map_limits_toward_open_space(
  labelled_map_base_limits
)
labelled_map_clusters <- find_map_inset_clusters(site_centres)

if (length(labelled_map_clusters) > 0L) {
  labelled_map_layout_plan <- plan_map_inset_layouts(
    site_centres,
    labelled_map_base_limits,
    labelled_map_clusters
  )
  labelled_map_layouts <- labelled_map_layout_plan$layouts
  labelled_map_limits <- labelled_map_layout_plan$map_limits
  labelled_map_zoom_limits_list <- lapply(
    seq_along(labelled_map_clusters),
    function(cluster_index) {
      cluster <- labelled_map_clusters[[cluster_index]]
      layout <- labelled_map_layouts[[cluster_index]]
      raw_limits <- list(xlim = cluster$xlim, ylim = cluster$ylim)
      fitted_limits <- fit_map_zoom_limits_to_inset(
        cluster$xlim,
        cluster$ylim,
        layout$position
      )
      blend_map_zoom_limits(raw_limits, fitted_limits)
    }
  )
} else {
  labelled_map_layouts <- list()
  labelled_map_limits <- labelled_map_base_limits
  labelled_map_zoom_limits_list <- list()
}

# Singular aliases retain compatibility with older code while the new list
# objects support zero, one, or two insets.
labelled_map_cluster <- if (length(labelled_map_clusters)) {
  labelled_map_clusters[[1]]
} else NULL
labelled_map_layout <- if (length(labelled_map_layouts)) {
  labelled_map_layouts[[1]]
} else NULL
labelled_map_zoom_limits <- if (length(labelled_map_zoom_limits_list)) {
  labelled_map_zoom_limits_list[[1]]
} else NULL
labelled_map_inset_sites <- unique(unlist(lapply(
  labelled_map_clusters,
  `[[`,
  "site_labels"
)))

labelled_map_subtitle <- paste0(
  scales::comma(sum(map_df$n_filtered, na.rm = TRUE)),
  " retained samples; ",
  scales::comma(n_sites_map),
  " site",
  ifelse(n_sites_map == 1, "", "s"),
  "; ",
  scales::comma(n_species_map),
  " species/group",
  ifelse(n_species_map == 1, "", "s")
)

species_legend_theme <- ggplot2::theme(
  legend.title = ggplot2::element_text(size = 11, face = "bold"),
  legend.text = ggplot2::element_text(size = 10),
  legend.key.width = grid::unit(0.78, "cm"),
  legend.key.height = grid::unit(0.78, "cm")
)
})

# ---- Map 10: species points with inset ----
rnr_map_chunks$map_10 <- quote({
#  Map 10: species/group points

if (is.null(labelled_map_cluster)) {
  species_labels_inset <- species_site_labels[0, , drop = FALSE]
} else {
  species_labels_inset <- species_site_labels %>%
    dplyr::filter(.data$site_label %in%
                    labelled_map_cluster$site_labels)
}

all_species_label_map <- make_map_base(
  labelled_map_limits,
  title = paste0(study_name, " — species/group locations"),
  subtitle = labelled_map_subtitle,
  caption = NULL,
  projection_crs = 4326
) +
  ggplot2::geom_point(
    data = species_site_labels,
    ggplot2::aes(
      x = .data$long,
      y = .data$lat,
      fill = .data$species_label,
      colour = .data$species_label
    ),
    shape = 21,
    size = map_labelled_point_size,
    stroke = 0.52,
    alpha = 1
  ) +
  ggplot2::scale_fill_manual(
    name = "Species/group",
    values = species_palette,
    breaks = species_levels,
    labels = species_legend_labels,
    drop = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = species_palette,
    breaks = species_levels,
    drop = FALSE,
    guide = "none"
  )

all_species_label_map <- all_species_label_map +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(
        shape = 21,
        size = map_species_legend_point_size,
        colour = unname(species_palette[species_levels])
      )
    )
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "left"
  ) +
  species_legend_theme

if (length(labelled_map_clusters) > 0L) for (
  cluster_index in seq_along(labelled_map_clusters)
) {
  current_cluster <- labelled_map_clusters[[cluster_index]]
  current_layout <- labelled_map_layouts[[cluster_index]]
  current_zoom_limits <- labelled_map_zoom_limits_list[[cluster_index]]
  species_labels_inset <- species_site_labels %>%
    dplyr::filter(.data$site_label %in% current_cluster$site_labels)

  species_zoom_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = aus_states,
      fill = "#F4F0E6",
      colour = "#6F7C82",
      linewidth = 0.28
    ) +
    ggplot2::geom_point(
      data = species_labels_inset,
      ggplot2::aes(
        x = .data$long,
        y = .data$lat,
        fill = .data$species_label,
        colour = .data$species_label
      ),
      shape = 21,
      size = map_labelled_point_size * map_inset_point_scale,
      stroke = 0.50,
      alpha = 1,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = species_palette, guide = "none") +
    ggplot2::scale_colour_manual(values = species_palette, guide = "none") +
    ggplot2::coord_sf(
      crs = sf::st_crs(4326),
      default_crs = sf::st_crs(4326),
      datum = sf::st_crs(4326),
      xlim = current_zoom_limits$xlim,
      ylim = current_zoom_limits$ylim,
      expand = FALSE,
      clip = "on"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        fill = "#EAF3F5",
        colour = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "#EAF3F5",
        colour = NA
      ),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  all_species_label_map <- add_labelled_map_inset(
    all_species_label_map,
    species_zoom_map,
    current_cluster,
    current_layout,
    current_zoom_limits
  )
}

all_species_label_dimensions <- make_labelled_map_dimensions(
  labelled_map_limits,
  n_labels = 0,
  minimum_width = 38
)

save_labelled_map_plot(
  "Map_10_All_species_points_with_inset",
  all_species_label_map,
  all_species_label_dimensions
)

#  Map 11: every site label -
})

# ---- Map 11: all site labels with inset ----
rnr_map_chunks$map_11 <- quote({
if (length(labelled_map_clusters) == 0L) {
  site_labels_main <- site_centres
  site_labels_inset <- site_centres[0, , drop = FALSE]
} else {
  site_labels_main <- site_centres %>%
    dplyr::filter(!.data$site_label %in%
                    labelled_map_inset_sites)
  site_labels_inset <- site_centres[0, , drop = FALSE]
}

all_site_label_map <- make_map_base(
  labelled_map_limits,
  title = paste0(study_name, " — all site labels"),
  subtitle = labelled_map_subtitle,
  caption = NULL,
  projection_crs = 4326
) +
  ggplot2::geom_point(
    data = site_centres,
    ggplot2::aes(
      x = .data$long,
      y = .data$lat,
      fill = .data$site_label,
      colour = .data$site_label
    ),
    shape = 21,
    size = map_labelled_point_size,
    stroke = 0.52,
    alpha = 1
  ) +
  ggplot2::scale_fill_manual(
    values = site_palette,
    breaks = site_levels,
    labels = site_legend_labels,
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::scale_colour_manual(
    values = site_palette,
    breaks = site_levels,
    drop = FALSE,
    guide = "none"
  ) +
  ggnewscale::new_scale_colour() +
  ggplot2::scale_colour_manual(
    values = site_label_palette,
    breaks = site_levels,
    drop = FALSE,
    guide = "none"
  )

all_site_label_map <- all_site_label_map +
  ggplot2::theme(legend.position = "none")

if (length(labelled_map_clusters) > 0L) for (
  cluster_index in seq_along(labelled_map_clusters)
) {
  current_cluster <- labelled_map_clusters[[cluster_index]]
  current_layout <- labelled_map_layouts[[cluster_index]]
  current_zoom_limits <- labelled_map_zoom_limits_list[[cluster_index]]
  site_labels_inset <- site_centres %>%
    dplyr::filter(.data$site_label %in% current_cluster$site_labels)

  site_zoom_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = aus_states,
      fill = "#F4F0E6",
      colour = "#6F7C82",
      linewidth = 0.28
    ) +
    ggplot2::geom_point(
      data = site_labels_inset,
      ggplot2::aes(
        x = .data$long,
        y = .data$lat,
        fill = .data$site_label,
        colour = .data$site_label
      ),
      shape = 21,
      size = map_labelled_point_size * map_inset_point_scale,
      stroke = 0.50,
      alpha = 1,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = site_palette,
      guide = "none"
    ) +
    ggplot2::scale_colour_manual(
      values = site_palette,
      guide = "none"
    ) +
    ggnewscale::new_scale_colour() +
    ggrepel::geom_label_repel(
      data = site_labels_inset,
      ggplot2::aes(
        x = .data$long,
        y = .data$lat,
        label = .data$display_site_count,
        colour = .data$site_label
      ),
      seed = map_seed,
      size = map_inset_label_text_size,
      fontface = "bold",
      fill = "white",
      alpha = 1,
      label.padding = grid::unit(0.10, "lines"),
      box.padding = grid::unit(0.22, "lines"),
      point.padding = grid::unit(0.08, "lines"),
      force = 1.0,
      force_pull = 3.5,
      segment.linewidth = 0.24,
      min.segment.length = 0,
      max.time = 25,
      max.iter = 100000,
      max.overlaps = Inf,
      xlim = current_zoom_limits$xlim,
      ylim = current_zoom_limits$ylim,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = site_label_palette,
      guide = "none"
    ) +
    ggplot2::coord_sf(
      crs = sf::st_crs(4326),
      default_crs = sf::st_crs(4326),
      datum = sf::st_crs(4326),
      xlim = current_zoom_limits$xlim,
      ylim = current_zoom_limits$ylim,
      expand = FALSE,
      clip = "on"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        fill = "#EAF3F5",
        colour = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "#EAF3F5",
        colour = NA
      ),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  all_site_label_map <- add_labelled_map_inset(
    all_site_label_map,
    site_zoom_map,
    current_cluster,
    current_layout,
    current_zoom_limits
  )
}

# Add direct labels last, so their opaque label boxes always sit above inset
# source frames and connector lines rather than being crossed by them.
all_site_label_map <- all_site_label_map +
  make_close_site_label_layer(
    site_labels_main,
    "display_site_count",
    limits = labelled_map_limits,
    text_size = max(3.6, label_text_size * 1.25),
    colour_field = "site_label"
  )

all_site_label_dimensions <- make_labelled_map_dimensions(
  labelled_map_limits,
  n_labels = n_sites_map,
  minimum_width = 38
)

save_labelled_map_plot(
  "Map_11_All_site_labels_with_inset",
  all_site_label_map,
  all_site_label_dimensions
)

#  Map 12: site labels with points coloured by species/group
})

# ---- Map 12: site labels coloured by species with inset ----
rnr_map_chunks$map_12 <- quote({
if (length(labelled_map_clusters) == 0L) {
  site_species_labels_main <- species_site_labels
  site_species_labels_inset <- species_site_labels[0, , drop = FALSE]
} else {
  site_species_labels_main <- species_site_labels %>%
    dplyr::filter(!.data$site_label %in%
                    labelled_map_inset_sites)
  site_species_labels_inset <- species_site_labels[0, , drop = FALSE]
}

site_by_species_map <- make_map_base(
  labelled_map_limits,
  title = paste0(study_name, " — site labels by species/group"),
  subtitle = labelled_map_subtitle,
  caption = NULL,
  projection_crs = 4326
) +
  ggplot2::geom_point(
    data = species_site_labels,
    ggplot2::aes(
      x = .data$long,
      y = .data$lat,
      fill = .data$species_label,
      colour = .data$species_label
    ),
    shape = 21,
    size = map_labelled_point_size,
    stroke = 0.52,
    alpha = 1
  ) +
  ggplot2::scale_fill_manual(
    name = "Species/group",
    values = species_palette,
    breaks = species_levels,
    labels = species_legend_labels,
    drop = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = species_palette,
    breaks = species_levels,
    drop = FALSE,
    guide = "none"
  ) +
  ggnewscale::new_scale_colour() +
  ggplot2::scale_colour_manual(
    values = species_label_palette,
    breaks = species_levels,
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(
        shape = 21,
        size = map_species_legend_point_size,
        colour = unname(species_palette[species_levels])
      )
    )
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "left"
  ) +
  species_legend_theme

if (length(labelled_map_clusters) > 0L) for (
  cluster_index in seq_along(labelled_map_clusters)
) {
  current_cluster <- labelled_map_clusters[[cluster_index]]
  current_layout <- labelled_map_layouts[[cluster_index]]
  current_zoom_limits <- labelled_map_zoom_limits_list[[cluster_index]]
  site_species_labels_inset <- species_site_labels %>%
    dplyr::filter(.data$site_label %in% current_cluster$site_labels)

  site_by_species_zoom_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = aus_states,
      fill = "#F4F0E6",
      colour = "#6F7C82",
      linewidth = 0.28
    ) +
    ggplot2::geom_point(
      data = site_species_labels_inset,
      ggplot2::aes(
        x = .data$long,
        y = .data$lat,
        fill = .data$species_label,
        colour = .data$species_label
      ),
      shape = 21,
      size = map_labelled_point_size * map_inset_point_scale,
      stroke = 0.50,
      alpha = 1,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = species_palette,
      guide = "none"
    ) +
    ggplot2::scale_colour_manual(
      values = species_palette,
      guide = "none"
    ) +
    ggnewscale::new_scale_colour() +
    ggrepel::geom_label_repel(
      data = site_species_labels_inset,
      ggplot2::aes(
        x = .data$long,
        y = .data$lat,
        label = .data$display_site_count,
        colour = .data$species_label
      ),
      seed = map_seed,
      size = map_inset_label_text_size,
      fontface = "bold",
      fill = "white",
      alpha = 1,
      label.padding = grid::unit(0.10, "lines"),
      box.padding = grid::unit(0.22, "lines"),
      point.padding = grid::unit(0.08, "lines"),
      force = 1.0,
      force_pull = 3.5,
      segment.linewidth = 0.24,
      min.segment.length = 0,
      max.time = 25,
      max.iter = 100000,
      max.overlaps = Inf,
      xlim = current_zoom_limits$xlim,
      ylim = current_zoom_limits$ylim,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = species_label_palette,
      guide = "none"
    ) +
    ggplot2::coord_sf(
      crs = sf::st_crs(4326),
      default_crs = sf::st_crs(4326),
      datum = sf::st_crs(4326),
      xlim = current_zoom_limits$xlim,
      ylim = current_zoom_limits$ylim,
      expand = FALSE,
      clip = "on"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        fill = "#EAF3F5",
        colour = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "#EAF3F5",
        colour = NA
      ),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  site_by_species_map <- add_labelled_map_inset(
    site_by_species_map,
    site_by_species_zoom_map,
    current_cluster,
    current_layout,
    current_zoom_limits
  )
}

# As above, draw labels after the inset components to preserve legibility where
# an inset source boundary lies close to a non-inset site.
site_by_species_map <- site_by_species_map +
  make_close_site_label_layer(
    site_species_labels_main,
    "display_site_count",
    limits = labelled_map_limits,
    text_size = max(3.6, label_text_size * 1.25),
    colour_field = "species_label"
  )

site_by_species_dimensions <- make_labelled_map_dimensions(
  labelled_map_limits,
  n_labels = n_populations_map,
  minimum_width = 38
)

save_labelled_map_plot(
  "Map_12_All_site_labels_by_species_with_inset",
  site_by_species_map,
  site_by_species_dimensions
)
})

# ---- Shared numbered-site key (Maps 13-14) ----
rnr_map_chunks$numbered <- quote({
#  Maps 13–14: numbered-site alternatives for dense labelling
# Keys use one column for <=20 sites and two balanced columns for larger
# datasets, while retaining north-to-south numbering.
site_number_key <- site_centres %>%
  dplyr::arrange(dplyr::desc(.data$lat), .data$long) %>%
  dplyr::mutate(
    site_number = dplyr::row_number(),
    site_number_label = as.character(.data$site_number)
  )

primary_species_by_site <- species_site_labels %>%
  dplyr::group_by(.data$site_label) %>%
  dplyr::arrange(dplyr::desc(.data$n_filtered), .data$species_label, .by_group = TRUE) %>%
  dplyr::slice_head(n = 1L) %>%
  dplyr::ungroup() %>%
  dplyr::select(.data$site_label, primary_species_label = .data$species_label)

site_number_key <- site_number_key %>%
  dplyr::left_join(primary_species_by_site, by = "site_label")

site_numbered_centres <- site_centres %>%
  dplyr::left_join(
    site_number_key %>%
      dplyr::select(
        .data$site_label,
        .data$site_number,
        .data$site_number_label,
        .data$primary_species_label
      ),
    by = "site_label"
  )

site_numbered_species_labels <- species_site_labels %>%
  dplyr::left_join(
    site_number_key %>%
      dplyr::select(.data$site_label, .data$site_number, .data$site_number_label),
    by = "site_label"
  )

numbered_map_dimensions <- make_labelled_map_dimensions(
  labelled_map_limits,
  n_labels = n_sites_map,
  minimum_width = 46
)
# Maps 13-14 do not need the full 3:2 landscape width. A 1.35:1 device uses
# the available space efficiently while leaving a sufficiently wide key panel
# for the enlarged two-column text.
numbered_map_dimensions["width"] <- max(
  46,
  min(82, unname(numbered_map_dimensions["height"]) * 1.35)
)
numbered_map_dimensions <- constrain_map_output_dimensions(
  numbered_map_dimensions
)
})

# ---- Map 13: site numbers with inset ----
rnr_map_chunks$map_13 <- quote({
numbered_site_key_plot <- make_numbered_site_key(
  site_number_key,
  colour_field = "site_label",
  palette = site_palette,
  title = "Site key (north to south)"
)

numbered_site_map <- make_numbered_site_map(
  point_data = site_numbered_centres,
  label_data = site_numbered_centres,
  colour_field = "site_label",
  palette = site_palette,
  title = paste0(study_name, " — site numbers")
)

save_labelled_map_plot(
  "Map_13_All_site_numbers_with_inset",
  numbered_site_map + numbered_site_key_plot +
    patchwork::plot_layout(widths = c(5.2, 2.8)),
  numbered_map_dimensions
)
})

# ---- Map 14: site numbers by species with inset ----
rnr_map_chunks$map_14 <- quote({
numbered_species_key_plot <- make_numbered_site_key(
  site_number_key,
  colour_field = "primary_species_label",
  palette = species_palette,
  title = "Site key (north to south)"
)

numbered_species_map <- make_numbered_site_map(
  point_data = site_numbered_species_labels,
  label_data = site_numbered_centres,
  colour_field = "species_label",
  palette = species_palette,
  title = paste0(study_name, " — site numbers by species/group"),
  label_colour_field = "primary_species_label"
)

save_labelled_map_plot(
  "Map_14_All_site_numbers_by_species_with_inset",
  numbered_species_map + numbered_species_key_plot +
    patchwork::plot_layout(widths = c(5.2, 2.8)),
  numbered_map_dimensions
)
})


#### 4. Map catalogue ####
rnr_map_catalogue <- data.frame(
  number = 0:14,
  key = c("context", "site_local", "site_ordered", "species_local", "species_ordered",
          "species_no_labels", "sampling_effort", "species_facets", "longitude_profile",
          "latitude_profile", "species_points", "site", "species", "site_numbers",
          "species_numbers"),
  stem = c("Map_00_Australia_context", "Map_01_Sites_local_labels", "Map_02_Sites_ordered_labels",
           "Map_03_Species_local_labels", "Map_04_Species_ordered_labels", "Map_05_Species_no_labels",
           "Map_06_Sampling_effort", "Map_07_Species_facets", "Map_08_Longitude_profile",
           "Map_09_Latitude_profile", "Map_10_All_species_points_with_inset",
           "Map_11_All_site_labels_with_inset", "Map_12_All_site_labels_by_species_with_inset",
           "Map_13_All_site_numbers_with_inset", "Map_14_All_site_numbers_by_species_with_inset"),
  needs = c("", "local_labels", "local_labels", "local_labels", "local_labels", "",
            "local_labels", "local_labels", "profile", "profile", "inset", "inset", "inset",
            "inset,numbered", "inset,numbered"),
  stringsAsFactors = FALSE
)
rnr_map_default_maps <- c("site", "species")

# Turns a maps argument into catalogue row numbers. Accepts keys ("site"),
# numbers (11 or "11"), stems ("Map_11_All_site_labels_with_inset" or "Map_11"),
# "all", "default" and "none".
resolve_rnr_maps <- function(maps) {
  if (is.null(maps) || isFALSE(maps) || length(maps) == 0L) return(integer())
  if (isTRUE(maps)) maps <- rnr_map_default_maps
  maps <- trimws(tolower(as.character(maps)))
  maps <- maps[!is.na(maps) & nzchar(maps)]
  if (any(maps %in% c("none", "false"))) return(integer())
  if (any(maps == "all")) return(seq_len(nrow(rnr_map_catalogue)))
  maps <- unlist(lapply(maps, function(x) if (x == "default") rnr_map_default_maps else x))
  rows <- vapply(maps, function(x) {
    hit <- which(tolower(rnr_map_catalogue$key) == x)
    if (!length(hit)) hit <- which(tolower(rnr_map_catalogue$stem) == x)
    if (!length(hit)) {
      number <- suppressWarnings(as.integer(sub("^map_?0*([0-9]+).*$", "\\1", x)))
      if (is.na(number)) number <- suppressWarnings(as.integer(x))
      if (!is.na(number)) hit <- which(rnr_map_catalogue$number == number)
    }
    if (length(hit)) hit[1] else NA_integer_
  }, integer(1))
  if (anyNA(rows)) {
    stop("Unknown map name(s): ", paste(maps[is.na(rows)], collapse = ", "),
         ". Use one or more of: ", paste(rnr_map_catalogue$key, collapse = ", "),
         ", a map number (0-14), \"all\" or \"none\".", call. = FALSE)
  }
  sort(unique(rows))
}

# Lists the available maps.
list_rnr_maps <- function() {
  out <- rnr_map_catalogue[, c("number", "key", "stem")]
  out$default <- out$key %in% rnr_map_default_maps
  out
}


#### 5. Main functions ####
# Builds base_map, divxlims and divylims (used by the scatterpie maps later in
# the workflow). Quick; run it even when no maps are drawn.
make_rnr_base_map <- function(envir = rnr_maps_env) {
  eval(rnr_map_chunks$base_map, envir = envir)
  invisible(get("base_map", envir = envir))
}

# Draws the requested maps into <output_dir>/plots/Maps.
#   maps          map names/numbers (see list_rnr_maps()); "all"; "none"
#   save_pdf      TRUE also saves PDFs (NULL keeps map_save_pdf)
#   dpi           PNG resolution for Maps 00-09 (NULL keeps map_dpi)
#   stop_on_error FALSE (default) reports a failing map and carries on with the
#                 other maps and the rest of the workflow; TRUE stops.
# Returns the saved file paths (invisibly).
run_rnr_maps <- function(maps = rnr_map_default_maps, save_pdf = NULL, dpi = NULL,
                         stop_on_error = FALSE, envir = rnr_maps_env) {
  rows <- resolve_rnr_maps(maps)
  if (!length(rows)) {
    message("MAPS: no maps requested; skipped.")
    return(invisible(character()))
  }
  if (!is.null(save_pdf)) assign("map_save_pdf", isTRUE(save_pdf), envir = envir)
  if (!is.null(dpi)) assign("map_dpi", dpi, envir = envir)

  selected <- rnr_map_catalogue[rows, , drop = FALSE]
  needs <- unique(unlist(strsplit(selected$needs[nzchar(selected$needs)], ",", fixed = TRUE)))
  message("MAPS: drawing ", nrow(selected), " map(s): ", paste(selected$stem, collapse = ", "))

  run_chunk <- function(chunk_name, label, fatal = FALSE) {
    started <- Sys.time()
    ok <- tryCatch({
      eval(rnr_map_chunks[[chunk_name]], envir = envir)
      TRUE
    }, error = function(e) {
      msg <- paste0("MAPS: ", label, " failed: ", conditionMessage(e))
      if (isTRUE(stop_on_error)) stop(msg, call. = FALSE)
      warning(msg, call. = FALSE, immediate. = TRUE)
      FALSE
    })
    if (ok) message("MAPS: ", label, " done (",
                    round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1), " s)")
    ok
  }

  # Maps 13-14 combine a map and a key with patchwork's "+" method, which
  # needs the patchwork namespace loaded first.
  if (any(selected$number %in% 13:14) && !requireNamespace("patchwork", quietly = TRUE)) {
    warning("MAPS: the patchwork package is needed for Maps 13-14; install it with install.packages('patchwork').",
            call. = FALSE)
  }
  if (!run_chunk("prep", "preparation")) {
    warning("MAPS: preparation failed, so no maps were drawn.", call. = FALSE)
    return(invisible(character()))
  }
  shared_ok <- c(local_labels = TRUE, profile = TRUE, inset = TRUE, numbered = TRUE)
  for (step in c("local_labels", "profile", "inset", "numbered")) {
    if (step %in% needs) shared_ok[[step]] <- run_chunk(step, paste0("shared step '", step, "'"))
  }
  for (i in seq_len(nrow(selected))) {
    map_needs <- strsplit(selected$needs[i], ",", fixed = TRUE)[[1]]
    map_needs <- map_needs[nzchar(map_needs)]
    if (length(map_needs) && !all(shared_ok[map_needs])) {
      warning("MAPS: ", selected$stem[i], " skipped because a shared step failed.", call. = FALSE)
      next
    }
    run_chunk(sprintf("map_%02d", selected$number[i]), selected$stem[i])
  }

  saved <- get0("saved_map_files", envir = envir, ifnotfound = character())
  maps_folder <- get0("maps_dir", envir = envir, ifnotfound = NA_character_)
  message("Map section complete: ", length(saved), " figure file(s) saved in ",
          normalizePath(maps_folder, winslash = "/", mustWork = FALSE))
  invisible(saved)
}
