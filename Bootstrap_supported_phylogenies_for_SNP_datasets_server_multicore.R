#!/usr/bin/env Rscript

### STAND-ALONE LINUX-SERVER PHYLOGENY WORKFLOW ###########################################
# This script is designed to run in a new R environment on the Linux server. It:
#   1. uses /home/richard/rrspecies/Astrotricha_all as the project directory;
#   2. installs missing R packages into a project-specific R_library folder;
#   3. reads species, site_col_name and species_col_name from the existing *_setup_variables.xlsx workbook;
#   4. reads the existing dms object from <species>/outputs_<site>_<species>/r_files/dms.RData;
#   5. recreates the site and species colour/shape palettes; and
#   6. runs the complete supported NJ/ML phylogeny workflow and writes the same SupportedPhylogenies output structure.
#
# Run from a terminal with:
#   cd /home/richard/rrspecies/Astrotricha_all
#   nohup Rscript --vanilla Bootstrap_supported_phylogenies_for_SNP_datasets_server_multicore.R > Bootstrap_supported_phylogenies_for_SNP_datasets_server_multicore.log 2>&1 &
#
# Monitor progress:
#       echo $! > Bootstrap_supported_phylogenies_for_SNP_datasets_server_multicore.pid
# 
#       tail -f Bootstrap_supported_phylogenies_for_SNP_datasets_server_multicore.log
# 
# To monitor the processes in another terminal:
#       ps -u richard -o pid,ppid,%cpu,%mem,etime,cmd --sort=-%cpu | head -20



#### A. SERVER PATH AND PACKAGE INSTALLATION ##############################################

project_dir <- "/home/richard/rrspecies/Astrotricha_all"


#main parameters to change:
# Bootstrap settings
phylo_nj_bootstrap_reps <- 1000L            # How many reps for NJ tree. Use 10-20 for testing and ~1000 or more for final analyses.
phylo_ml_bootstrap_reps <- 1000L            # How many reps for ML tree. Use 10-20 for testing and ~1000 or more for final analyses.
phylo_bootstrap_seed <- 123L
phylo_support_to_plot <- 50                 # only show values of support on figures that are above this value 

# Linux parallel processing
phylo_use_parallel <- FALSE                  # Set to TRUE to use multiple cores (for server runs)
phylo_requested_cores <- 4L                 # Set this to the number of CPU cores allocated to this analysis.
phylo_parallel_bootstraps <- TRUE           # Runs independent NJ and ML bootstrap replicates concurrently.
phylo_parallel_edge_starts <- TRUE          # Runs the five Stage 1 and Stage 3 branch-length starting values concurrently.





setup_file <- NULL                         # NULL = find the single *_setup_variables.xlsx file in project_dir.
dms_file <- NULL                           # NULL = construct the expected dms path from the setup workbook.
install_missing_packages <- TRUE
install_optional_ragg <- FALSE             # ragg improves PNG rendering but may need Linux system libraries; the workflow works without it.
cran_repository <- "https://cloud.r-project.org"
local_R_library <- file.path(project_dir, "R_library")

if (!dir.exists(project_dir)) stop("Project directory does not exist: ", project_dir, call.=FALSE)
dir.create(local_R_library, recursive=TRUE, showWarnings=FALSE)
.libPaths(unique(c(local_R_library, .libPaths())))
options(repos=c(CRAN=cran_repository), timeout=max(600, getOption("timeout", 60)))
setwd(project_dir)

required_packages <- c("openxlsx", "ape", "phangorn", "phytools", "colorspace")
optional_packages <- if (isTRUE(install_optional_ragg)) "ragg" else character(0)
missing_packages <- c(required_packages, optional_packages)[!vapply(c(required_packages, optional_packages), requireNamespace, logical(1), quietly=TRUE)]

if (length(missing_packages)>0L && isTRUE(install_missing_packages)) {
  install_cores <- suppressWarnings(parallel::detectCores(logical=FALSE)); if (!is.finite(install_cores)) install_cores <- 1L
  install_cores <- max(1L, as.integer(install_cores)-1L)
  message("[Package setup] Installing missing packages into: ", local_R_library)
  message("[Package setup] Packages: ", paste(missing_packages, collapse=", "))
  install.packages(missing_packages, lib=local_R_library, dependencies=c("Depends", "Imports", "LinkingTo"), Ncpus=install_cores)
}

still_missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(still_missing)>0L) stop("Required packages are still unavailable after installation: ", paste(still_missing, collapse=", "),
                                  "\nCheck the installation messages. The server may require Linux build tools or system libraries.", call.=FALSE)
message("[Package setup] Required R packages are available. Library paths: ", paste(.libPaths(), collapse=" | "))


#### B. READ SETUP VARIABLES AND THE EXISTING DMS ########################################

if (is.null(setup_file)) {
  expected_setup <- file.path(project_dir, paste0(basename(project_dir), "_setup_variables.xlsx"))
  setup_candidates <- list.files(project_dir, pattern="_setup_variables\\.xlsx$", full.names=TRUE, ignore.case=TRUE)
  if (file.exists(expected_setup)) setup_file <- expected_setup
  else if (length(setup_candidates)==1L) setup_file <- setup_candidates
  else if (length(setup_candidates)==0L) stop("No *_setup_variables.xlsx workbook was found in: ", project_dir, call.=FALSE)
  else stop("More than one setup workbook was found. Set setup_file manually at the top of this script:\n", paste(setup_candidates, collapse="\n"), call.=FALSE)
}
if (!file.exists(setup_file)) stop("Setup workbook not found: ", setup_file, call.=FALSE)

setup_variables <- openxlsx::read.xlsx(setup_file, colNames=TRUE)
read_setup_value <- function(row_number, description) {
  if (nrow(setup_variables)<row_number || ncol(setup_variables)<2L) stop("The setup workbook does not contain the expected row for ", description, ".", call.=FALSE)
  value <- trimws(as.character(setup_variables[row_number, 2]))
  if (length(value)!=1L || is.na(value) || value=="") stop("The setup workbook value is blank for ", description, " (row ", row_number, ", column 2).", call.=FALSE)
  value
}

species <- read_setup_value(2L, "species")
species_col_name <- read_setup_value(5L, "species_col_name")
site_col_name <- read_setup_value(6L, "site_col_name")
output_root <- file.path(project_dir, species, paste0("outputs_", site_col_name, "_", species_col_name))

if (is.null(dms_file)) dms_file <- file.path(output_root, "r_files", "dms.RData")
if (!file.exists(dms_file)) {
  dms_candidates <- list.files(project_dir, pattern="^dms\\.(RData|rds)$", recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
  if (length(dms_candidates)==1L) { dms_file <- dms_candidates; warning("Expected dms path was not found; using the only dms file located: ", dms_file, call.=FALSE) }
  else stop("The expected dms file was not found: ", dms_file, "\nLocated alternatives:\n", if (length(dms_candidates)) paste(dms_candidates, collapse="\n") else "None", call.=FALSE)
}

message("[DMS setup] Reading: ", dms_file)
dms <- readRDS(dms_file)
if (!is.list(dms) || is.null(dms$gt) || is.null(dms$meta$analyses)) stop("The loaded object is not a valid dms with dms$gt and dms$meta$analyses.", call.=FALSE)

# dms$meta$analyses may be stored as either a matrix or a data frame. Convert it once
# so named metadata columns behave consistently throughout the stand-alone workflow.
dms$meta$analyses <- as.data.frame(dms$meta$analyses, stringsAsFactors=FALSE, check.names=FALSE)
metadata_columns <- colnames(dms$meta$analyses)
missing_metadata_columns <- setdiff(c(site_col_name, species_col_name), metadata_columns)
if (length(missing_metadata_columns)>0L) stop("The configured site/species columns were not found in dms$meta$analyses: ",
                                                paste(missing_metadata_columns, collapse=", "),
                                                "\nAvailable metadata columns are: ", paste(metadata_columns, collapse=", "), call.=FALSE)
message("[DMS setup] Loaded ", nrow(dms$gt), " samples × ", format(ncol(dms$gt), big.mark=","), " loci.")
message("[DMS setup] Species column: ", species_col_name, " | site column: ", site_col_name, " | output root: ", output_root)


#### C. RECREATE SPECIES AND SITE COLOUR/SHAPE PALETTES ##################################

clean_palette_group <- function(x) { x <- as.character(x); x[is.na(x) | trimws(x)==""] <- "Missing"; x }
species_groups <- clean_palette_group(dms$meta$analyses[, species_col_name, drop=TRUE])
site_groups <- clean_palette_group(dms$meta$analyses[, site_col_name, drop=TRUE])

species_levels <- sort(unique(species_groups))
sp_colours <- stats::setNames(colorspace::rainbow_hcl(length(species_levels)), species_levels)
sp_shapes <- stats::setNames(rep(21:25, length.out=length(species_levels)), species_levels)

# Retain the original case-insensitive colour rules used in the main workflow.
species_names_lower <- tolower(gsub("_", " ", names(sp_colours), fixed=TRUE))
sp_colours[grepl("local", species_names_lower)] <- "blue"
sp_colours[grepl("future|external", species_names_lower)] <- "orange"
sp_colours[grepl("nursery", species_names_lower)] <- "black"
sp_colours[species_names_lower %in% c("in situ", "insitu")] <- grDevices::adjustcolor("black", alpha.f=0.10)
if ("Missing" %in% names(sp_colours)) sp_colours["Missing"] <- "grey50"

site_levels <- unique(site_groups)
site_levels <- c(setdiff(site_levels, c("no_geo_data", "Missing")), intersect(c("no_geo_data", "Missing"), site_levels))
regular_site_levels <- setdiff(site_levels, c("no_geo_data", "Missing"))
site_colours <- stats::setNames(character(length(site_levels)), site_levels)
if (length(regular_site_levels)>0L) site_colours[regular_site_levels] <- colorspace::rainbow_hcl(length(regular_site_levels))
if ("no_geo_data" %in% site_levels) site_colours["no_geo_data"] <- "grey30"
if ("Missing" %in% site_levels) site_colours["Missing"] <- "grey50"
site_shapes <- stats::setNames(rep(c(21:25,3,4,8,12), length.out=length(site_levels)), site_levels)
if ("Missing" %in% names(site_shapes)) site_shapes["Missing"] <- 4

palette_output_dir <- file.path(output_root, "plots", "SupportedPhylogenies")
dir.create(palette_output_dir, recursive=TRUE, showWarnings=FALSE)
utils::write.csv(data.frame(species=names(sp_colours), colour=unname(sp_colours), shape=unname(sp_shapes), stringsAsFactors=FALSE),
                 file.path(palette_output_dir, "Phylogenetic_species_palette.csv"), row.names=FALSE)
utils::write.csv(data.frame(site=names(site_colours), colour=unname(site_colours), shape=unname(site_shapes), stringsAsFactors=FALSE),
                 file.path(palette_output_dir, "Phylogenetic_site_palette.csv"), row.names=FALSE)
message("[Palette setup] Assigned ", length(sp_colours), " species colours/shapes and ", length(site_colours), " site colours/shapes.")


### BOOTSTRAP-SUPPORTED NJ AND ML PHYLOGENIES ############################################
# SUMMARY
# Initial workflow developed by R.Dimon July 2026 using GPT5.6Sol for SNP genotype data
# This workflow generates neighbour-joining (NJ) and maximum-likelihood (ML) trees from a filtered SNP genotype matrix. It calculates locus-bootstrap support,
# produces site- and species-labelled tree figures, compares the NJ and ML topologies using tanglegrams, and generates missing-data diagnostic plots.


# INPUT DATA
# The genotype matrix is taken from dms$gt, with samples in rows and SNP loci in columns. Genotypes are coded as:
#   0 = homozygous reference genotype
#   1 = heterozygous genotype
#   2 = homozygous alternate genotype
#   NA = missing genotype call
# Sample names and grouping information are taken from dms$meta$analyses. The script checks that the genotype matrix and metadata contain the same samples and
# that sample names are unique. Loci that are completely missing or have fewer than two observed genotype states (fixed) are removed before analysis (should already be filtered out).
# The same filtered set of polymorphic SNP loci is used for the NJ tree, ML tree, bootstrap analyses and missing-data summaries.


# NEIGHBOUR-JOINING ANALYSIS
# The NJ tree is calculated from pairwise Euclidean distances among samples. These distances summarise how different the 0/1/2 genotype profiles are across all analysed SNP loci.
# NJ branch lengths are measured in Euclidean genotype-distance units (i.e., do not represent time, generations or substitutions per site).
# Long branches of one sample may indicate strong genetic differentiation, an unusual sample, or a possible effect of missing genotype data.
# Note, as the SNP genotype matrix uses only polymorphic loci, the terminal benches may appear longer than what you may expect with other datasets (i.e., concatenated markers)


# NJ LOCUS BOOTSTRAPPING
# NJ bootstrap support is calculated by repeatedly resampling SNP loci with replacement. Each replicate contains the same number of loci as the original
# dataset, but some loci may be sampled more than once while others may not be sampled.
# All samples are retained in every replicate. A new Euclidean-distance matrix and NJ tree are calculated from each resampled dataset. Support for an
# internal split is the percentage of bootstrap trees that recover the same split as the reference NJ tree.


# MAXIMUM-LIKELIHOOD ANALYSIS
# The ML analysis treats genotypes 0, 1 and 2 as three possible character states. 
# Missing genotypes are coded as unknown, meaning they are allowed to match any of the three states during likelihood calculations.
# Ascertainment-bias correction (+ASC) is used because only polymorphic SNP loci are analysed. 
# This correction accounts for the fact that monomorphic loci were removed before the likelihood analysis.
# The NJ topology is used as the starting tree for the ML analysis because it provides a fast, data-based starting arrangement of the samples. 
# NJ branch lengths are discarded before ML optimisation because NJ and ML branch lengths are measured in different units.


# ML MODEL OPTIONS
# phylo_ml_model controls how changes among genotype states are modelled:
#
#   "ER"
#   Uses an equal-rates model. Changes between all pairs of states are allowed at the same rate: 0 <-> 1, 0 <-> 2 and 1 <-> 2.
#
#   "ORDERED"
#   Treats genotype states as a stepwise series: 0 <-> 1 <-> 2. Direct changes between states 0 and 2 are not allowed.
#
# ER and ORDERED analyses can be compared using their final log-likelihoods, tree topology, branch lengths and bootstrap support, provided they use the same samples and SNP loci.


# ML OPTIMISATION
# ML branch lengths and topology are selected by maximising the likelihood of the observed genotype data under the chosen model. A higher log-likelihood
# indicates a better fit between the model, branch lengths, topology and observed data.



# The analysis uses three optimisation stages. These stages separate the estimation of branch lengths from the search for a better tree topology and
# provide checks that the final result is not strongly influenced by the initial branch-length values.
#
# Stage 1: Initial fixed-topology branch-length optimisation
#   The topology of the NJ tree is retained as the starting arrangement of the samples, but the NJ branch lengths are discarded because they are measured
#   in Euclidean genotype-distance units rather than ML model units. The same starting length is initially assigned to every branch. The analysis is repeated using each value listed in phylo_ml_edge_start_values, which
#   by default is: 0.0001, 0.001, 0.01, 0.05 and 0.10
#   These values are only initial guesses. During each optimisation run, individual branches are free to increase or decrease independently. The tree topology remains fixed throughout Stage 1, so only branch lengths are adjusted.
#   The final log-likelihood from each starting value is compared, and the fit with the highest likelihood is retained as the starting point for the topology search. 
#   Using several starting values helps detect cases where branch-length optimisation may become trapped in a poorer local solution.
#
# Stage 2: ML topology search
#   The best fixed-topology fit from Stage 1 is used to begin the topology search. During this stage, the arrangement of branches is allowed to change and branch lengths are reoptimised as alternative topologies are examined.
#   The type of topology search is controlled by phylo_ml_rearrangement:
#       "NNI"        tests local nearest-neighbour rearrangements around internal branches;
#       "stochastic" uses broader random perturbations followed by likelihood optimisation;
#       "ratchet"    (current default) - temporarily changes the influence of subsets of loci to encourage movement between different parts of topology space.
#   For each candidate arrangement, the likelihood of the observed genotype-state data is recalculated under the selected ML model and +ASC correction.
#   Topological changes that improve the likelihood are retained, and the search continues until the selected search method no longer finds a better fit.
#   The result of Stage 2 is therefore the highest-likelihood combination of topology and branch lengths found during that search. It is not guaranteed to
#   be the absolute global optimum, particularly for large datasets, but stochastic and ratchet searches generally explore more alternatives than a single local NNI search.
#
# Stage 3: Final fixed-topology branch-length refinement
#   The best topology found in Stage 2 is then held fixed. Its branch lengths are discarded and re-estimated again from each value in phylo_ml_edge_start_values.
#   As in Stage 1, every branch initially receives the same starting value, but branches are then optimised independently. No further topological changes
#   are allowed during this stage. This provides a final check that the branch lengths associated with the selected topology are stable and are not overly
#   dependent on the branch lengths inherited from the topology-search process.
#   The highest-likelihood Stage 3 fit is compared with the original Stage 2 topology-search fit. The refined fit is retained only when its likelihood is
#   equal to or better than the Stage 2 fit within phylo_ml_loglik_tolerance.
#   If the Stage 3 refinement produces a lower likelihood beyond this tolerance, the script keeps the Stage 2 topology-search fit instead and reports a
#   warning. This prevents a poorer fixed-topology refinement from replacing a better result already found during the topology search.
#
#   Results from all starting values are written to ML_staged_optimisation_diagnostics.csv. Similar final likelihoods and branch-length summaries among
#   the alternative starts suggest stable optimisation. Large differences among starts may indicate sensitivity to starting conditions or convergence on different local optima.


# ML LOCUS BOOTSTRAPPING
# ML bootstrap support is calculated by repeatedly resampling SNP loci with replacement. A new ML tree is fitted for each replicate using the selected model
# and +ASC correction. Branch lengths are optimised and NNI rearrangements are used to improve each replicate topology.
# The replicate trees are used to calculate support for the final reference ML tree (They do not replace the final reference topology or its branch lengths).
# ML bootstrapping is much slower than NJ bootstrapping because a likelihood analysis and topology improvement are carried out for every replicate.


# INTERPRETATION OF BOOTSTRAP SUPPORT
# Bootstrap support describes how consistently an internal split is recovered when SNP loci are resampled. For example, support of 95 means that the split
# was recovered in 95% of the bootstrap replicate trees.
# Bootstrap support measures resampling stability. It is not the probability that a clade is biologically correct.
# Approximate descriptive categories are:
#   >= 95%   very strong support
#   80-94%   strong support
#   70-79%   moderate support
#   50-69%   weak support
#   < 50%    poorly supported
# phylo_support_to_plot controls which support values are printed on the figures. Values below this threshold are still calculated and retained in the tree objects and support tables.
# Small numbers of bootstrap replicates are suitable only for testing the script. With 10 replicates, support values can change only in steps of 10%.
# Larger numbers of replicates are required for final analyses, and support stability should be checked as the number of replicates increases.


# ROOTING
# NJ and ML analyses produce unrooted trees. Rooting is applied only to copies used for plotting, while the saved analytical tree objects remain unrooted.
# One or more outgroup taxa can be entered in phylo_outgroup_taxa using exact values from species_col_name. All samples assigned to those taxa are used as
# outgroup tips. When several outgroup samples are used, they must form a single monophyletic group.
# When no valid outgroup is supplied, the plotting copy is midpoint rooted. Midpoint rooting places the root halfway along the longest path between two tips.
# It is used only to improve visual presentation and does not demonstrate the true evolutionary direction of the tree.


# STAND-ALONE NJ AND ML FIGURES
# Separate rectangular NJ and ML trees are produced using site and species groupings.
# Site-coloured figures use tip labels in the following format: species | site | sample accession
# Species-coloured figures use: species | sample accession
# Expanded labels are added only to plotting copies. The analytical tree objects, Newick files and bootstrap trees retain the original accession-only labels.
# phylo_tree_show_branch_lengths = TRUE produces phylograms using the estimated branch lengths. Setting it to FALSE produces cladogram-style figures with
# even visual branch spacing while retaining the original branch lengths in the saved tree objects.
# Figure width, height, label size and legend layout are adjusted according to the number of samples and the length of the tip labels. Trees are saved as vector PDFs and low-res PNG files


# TANGLEGRAMS
# A tanglegram places the NJ tree on the left and the ML tree on the right. Matching samples are connected by lines, allowing differences in sample placement between the two topologies to be seen directly.
# Internal nodes are rotated without changing tree topology to reduce the number of crossing links. Several starting arrangements can be tested, and the
# arrangement with the lowest tip-order mismatch score is retained.
# If automatic rotation fails, the script uses the original unrotated arrangement rather than stopping the analysis. Rotation scores and any errors are written to diagnostic output files.
# phylo_tanglegram_use_branch_lengths = FALSE gives all branches equal visual lengths and is generally easier for comparing topology. Setting it to TRUE
# retains the original NJ and ML branch-length patterns, although the two scales must not be compared numerically.
# Site-coloured tanglegrams use species | site | sample accession labels. Species-coloured tanglegrams use species | sample accession labels.


# MISSING-DATA DIAGNOSTICS
# Missing data are calculated separately for each sample as: 100 x number of missing genotype calls / number of analysed polymorphic loci
# Missing-data NJ, ML and tanglegram figures colour the labels, symbols and connecting lines according to the percentage of missing genotype calls. The missing-data % is also added to each plotting label.
# These figures help determine whether unusual tree placement or long terminal branches are associated with poor genotype completeness. 
# Missingness is used only as a plotting diagnostic and does not alter the NJ or ML tree calculations.


# CHECKPOINTS AND REPLOTTING
# NJ and ML bootstrap checkpoints are saved at regular intervals. If the analysis is interrupted, it can continue from the latest compatible checkpoint.
# phylo_reset_nj_checkpoint and phylo_reset_ml_checkpoint can be used to remove an old or incompatible checkpoint before restarting an analysis.
# phylo_replot_only = TRUE reloads previously saved NJ and ML trees, bootstrap trees and the ML model fit. This allows the figures and summary tables to be
# regenerated without repeating the tree estimation or bootstrap analyses.


# OUTPUTS
# The workflow produces:
#   - site- and species-coloured NJ and ML trees in PDF and PNG format;
#   - site-, species- and missing-data tanglegrams;
#   - supported NJ and ML trees in Newick and RDS formats;
#   - optional NJ and ML bootstrap replicate trees;
#   - the final ML model fit;
#   - ML optimisation diagnostics;
#   - bootstrap-support tables;
#   - tanglegram rotation diagnostics;
#   - sample missing-data summaries;
#   - tip metadata, analysis summaries and session information.


# REQUIRED PACKAGES
# openxlsx is used to read the setup workbook. ape, phangorn, phytools and colorspace are required; ragg is optional for PNG output.


# REQUIRED OBJECTS FROM THE MAIN WORKFLOW
# dms$gt               = filtered genotype matrix, with samples in rows and loci in columns
# dms$meta$analyses    = sample metadata
# species              = analysis name used in figure titles and output paths
# site_col_name        = name of the metadata column containing sampling-site groups
# species_col_name     = name of the metadata column containing species or taxon groups
# site_colours         = named colours for sampling sites
# site_shapes          = named plotting symbols for sampling sites
# sp_colours           = named colours for species or taxa
# sp_shapes            = named plotting symbols for species or taxa


#### 1. PHYLOGENY USER SETTINGS #########################################################
# Defaults below run NJ, ML and tanglegrams. The current bootstrap settings use 1000 NJ replicates and 10 ML replicates.
# For a quick server test, temporarily use 3-10 replicates and consider phylo_ml_rearrangement <- "NNI".

run_supported_phylogenies <- TRUE
phylo_replot_only <- FALSE                 # TRUE reloads saved supported trees/bootstrap trees and redraws outputs without rerunning NJ/ML estimation.

# Analyses and outputs
phylo_run_nj <- TRUE
phylo_run_ml <- TRUE
phylo_run_tanglegrams <- TRUE
phylo_generate_missing_data_plots <- TRUE
phylo_save_bootstrap_trees <- TRUE
phylo_save_collapsed_species_trees <- TRUE

# Metadata/output overrides
phylo_sample_col <- NULL                   # NULL = auto-detect; otherwise exact metadata column name.
phylo_output_dir <- NULL                   # NULL = species/outputs_<site>_<species>/plots/SupportedPhylogenies.



# ML settings
phylo_ml_model <- "ER"                     # "ER" or "ORDERED".
phylo_ml_rearrangement <- "NNI"            # "NNI", "stochastic" or "ratchet". NNI avoids the ratchet finalisation error seen on this server.
phylo_ml_edge_start_values <- c(1e-4, 1e-3, 0.01, 0.05, 0.10)
phylo_ml_edge_epsilon <- 1e-10
phylo_ml_edge_maxit <- 100L
phylo_ml_loglik_tolerance <- 1e-6
phylo_ml_optim_trace <- 1L

# Rooting
phylo_outgroup_taxa <- character(0)        # Exact values from species_col_name; empty = midpoint rooting for display only.

# Stand-alone NJ and ML tree plotting
phylo_tree_colour_sample_labels <- TRUE
phylo_tree_show_branch_lengths <- TRUE
phylo_tree_tip_label_separator <- " | "
phylo_tree_legend_max_columns <- 6L
phylo_tree_height_per_tip_cm <- 0.50
phylo_tree_min_height_cm <- 70
phylo_tree_max_height_cm <- 220
phylo_tree_min_width_cm <- 85
phylo_tree_max_width_cm <- 150
phylo_tree_png_res <- 300L
phylo_tree_large_png_res <- 120L
phylo_tree_large_figure_threshold_cm <- 100

# Collapsed monophyletic-species figures
phylo_collapsed_triangle_length_fraction <- 0.04   # Triangle length as a proportion of the plotted tree width.
phylo_collapsed_triangle_half_height <- 0.5        # Triangle extends this far above and below its branch.
phylo_collapsed_triangle_border_width <- 1.1
phylo_collapsed_label_scale <- 1.15                 # Make labels in collapsed trees slightly larger.

# Tanglegram plotting and order optimisation
phylo_tanglegram_use_branch_lengths <- FALSE
phylo_tanglegram_optimisation_starts <- 5L
phylo_tanglegram_random_fraction <- 0.40
phylo_tanglegram_rotation_seed <- 456L
phylo_tanglegram_methods <- c("pre", "post")
phylo_tanglegram_rotate_multifurcations <- TRUE
phylo_tanglegram_link_width <- 1.8
phylo_tanglegram_link_alpha <- 0.75
phylo_tanglegram_tree_width <- 0.85
phylo_tanglegram_tip_line_width <- 0.45
phylo_tanglegram_tip_line_type <- "dotted"
phylo_tanglegram_height_per_tip_cm <- 0.48
phylo_tanglegram_min_height_cm <- 80
phylo_tanglegram_max_height_cm <- 220
phylo_tanglegram_min_width_cm <- 120
phylo_tanglegram_max_width_cm <- 200

# Missing-data figures
phylo_missing_palette_name <- "YlOrRd"
phylo_missing_high_values_dark <- TRUE
phylo_missing_label_digits <- 2L

# Progress/checkpoints
phylo_progress_every <- 5L
phylo_checkpoint_every <- 5L
phylo_resume_bootstraps <- TRUE
phylo_progress_recent_n <- 5L
phylo_reset_nj_checkpoint <- FALSE
phylo_reset_ml_checkpoint <- FALSE


#### 2. SETUP, VALIDATION AND GENERAL HELPERS #############################################

if (run_supported_phylogenies) {
  
  #### 2.1 Packages, settings and output directory
  
  required_packages <- c("ape", "phangorn", "phytools", "colorspace")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing_packages)>0L) stop("Install required packages: ", paste(missing_packages, collapse=", "), call.=FALSE)

  # Use forked worker processes on Linux for independent bootstrap replicates and ML edge-start fits.
  if (!is.numeric(phylo_requested_cores) || length(phylo_requested_cores)!=1L || !is.finite(phylo_requested_cores) || phylo_requested_cores<1L) stop("phylo_requested_cores must be a positive integer.", call.=FALSE)
  phylo_requested_cores <- as.integer(phylo_requested_cores)
  phylo_detected_cores <- suppressWarnings(parallel::detectCores(logical=FALSE))
  if (!is.finite(phylo_detected_cores) || phylo_detected_cores<1L) phylo_detected_cores <- 1L
  phylo_scheduler_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset=NA_character_)))
  phylo_available_cores <- if (is.finite(phylo_scheduler_cores) && phylo_scheduler_cores>0L) phylo_scheduler_cores else phylo_detected_cores
  phylo_n_cores <- if (isTRUE(phylo_use_parallel) && .Platform$OS.type!="windows") max(1L,min(phylo_requested_cores,phylo_available_cores)) else 1L

  # Prevent each forked worker from starting additional BLAS/OpenMP threads and oversubscribing the server.
  Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
  message(sprintf("[Parallel setup] Using %d worker core(s); requested=%d, detected physical cores=%d%s.", phylo_n_cores, phylo_requested_cores, phylo_detected_cores,
                  if (is.finite(phylo_scheduler_cores) && phylo_scheduler_cores>0L) paste0(", scheduler allocation=",phylo_scheduler_cores) else ""))
  if (!phylo_ml_model %in% c("ORDERED", "ER")) stop("phylo_ml_model must be 'ORDERED' or 'ER'.", call.=FALSE)
  if (!phylo_ml_rearrangement %in% c("NNI", "stochastic", "ratchet")) stop("phylo_ml_rearrangement must be 'NNI', 'stochastic' or 'ratchet'.", call.=FALSE)
  if (phylo_run_nj && phylo_nj_bootstrap_reps<1L) stop("phylo_nj_bootstrap_reps must be at least 1.", call.=FALSE)
  if (phylo_run_ml && phylo_ml_bootstrap_reps<1L) stop("phylo_ml_bootstrap_reps must be at least 1.", call.=FALSE)
  if (phylo_tanglegram_optimisation_starts<1L) stop("phylo_tanglegram_optimisation_starts must be at least 1.", call.=FALSE)
  if (phylo_progress_every<1L || phylo_checkpoint_every<1L) stop("Progress and checkpoint intervals must be at least 1.", call.=FALSE)
  if (phylo_support_to_plot<0 || phylo_support_to_plot>100) stop("phylo_support_to_plot must be between 0 and 100.", call.=FALSE)
  if (any(!is.finite(phylo_ml_edge_start_values)) || any(phylo_ml_edge_start_values<=0)) stop("phylo_ml_edge_start_values must contain positive finite values.", call.=FALSE)
  
  configured_outgroup_taxa <- unique(as.character(phylo_outgroup_taxa))
  configured_outgroup_taxa <- configured_outgroup_taxa[!is.na(configured_outgroup_taxa) & trimws(configured_outgroup_taxa)!=""]
  
  phylo_dir <- if (is.null(phylo_output_dir)) file.path(species, paste0("outputs_", site_col_name, "_", species_col_name), "plots", "SupportedPhylogenies") else as.character(phylo_output_dir)[1]
  dir.create(phylo_dir, recursive=TRUE, showWarnings=FALSE)
  
  
  #### 2.2 Genotype matrix and metadata
  
  phylo_gt <- as.matrix(dms$gt); storage.mode(phylo_gt) <- "numeric"
  phylo_meta <- as.data.frame(dms$meta$analyses, stringsAsFactors=FALSE)
  
  if (is.null(phylo_sample_col)) {
    sample_col_candidates <- c("sample", "sample_names", "Sample", "SampleID", "sample_id", "individual", "id", "ID")
    phylo_sample_col <- sample_col_candidates[sample_col_candidates %in% names(phylo_meta)][1]
    if (is.na(phylo_sample_col)) stop("No sample-name column was found in dms$meta$analyses. Set phylo_sample_col manually.", call.=FALSE)
  } else {
    phylo_sample_col <- as.character(phylo_sample_col)[1]
    if (!phylo_sample_col %in% names(phylo_meta)) stop("phylo_sample_col was not found in dms$meta$analyses: ", phylo_sample_col, call.=FALSE)
  }
  
  if (!all(c(site_col_name, species_col_name) %in% colnames(phylo_meta))) stop("The selected site or species metadata column is missing.", call.=FALSE)
  phylo_meta$.phylo_sample <- as.character(phylo_meta[[phylo_sample_col]])
  
  if (is.null(rownames(phylo_gt)) || anyNA(rownames(phylo_gt)) || any(rownames(phylo_gt)=="")) {
    if (!is.null(dms$sample_names) && length(dms$sample_names)==nrow(phylo_gt)) rownames(phylo_gt) <- as.character(dms$sample_names)
    else if (nrow(phylo_meta)==nrow(phylo_gt)) rownames(phylo_gt) <- phylo_meta$.phylo_sample
    else stop("Sample names could not be assigned to dms$gt.", call.=FALSE)
  }
  
  if (anyDuplicated(rownames(phylo_gt))) stop("Duplicate sample names were found in dms$gt.", call.=FALSE)
  if (anyDuplicated(phylo_meta$.phylo_sample)) stop("Duplicate sample names were found in dms$meta$analyses.", call.=FALSE)
  
  phylo_meta <- phylo_meta[match(rownames(phylo_gt), phylo_meta$.phylo_sample), , drop=FALSE]
  if (anyNA(phylo_meta$.phylo_sample)) stop("Samples missing from metadata: ", paste(head(rownames(phylo_gt)[is.na(phylo_meta$.phylo_sample)], 20), collapse=", "), call.=FALSE)
  if (!identical(rownames(phylo_gt), phylo_meta$.phylo_sample)) stop("Metadata could not be aligned to the genotype matrix.", call.=FALSE)
  if (nrow(phylo_gt)<4L) stop("At least four samples are required.", call.=FALSE)
  if (is.null(colnames(phylo_gt))) colnames(phylo_gt) <- paste0("Locus_", seq_len(ncol(phylo_gt)))
  if (any(!is.na(phylo_gt) & !phylo_gt %in% c(0,1,2))) stop("dms$gt contains values other than 0, 1, 2 or NA.", call.=FALSE)
  if (any(rowSums(!is.na(phylo_gt))==0L)) stop("One or more samples contain no non-missing genotype calls.", call.=FALSE)
  
  phylo_n_states <- apply(phylo_gt, 2, function(z) length(unique(z[!is.na(z)])))
  phylo_keep_loci <- phylo_n_states>=2L
  message(sprintf("[Phylogeny setup] Dataset ready: %d samples × %s polymorphic loci; %s non-variable/all-missing loci removed.", nrow(phylo_gt), format(sum(phylo_keep_loci), big.mark=","), format(sum(!phylo_keep_loci), big.mark=",")))
  phylo_gt <- phylo_gt[, phylo_keep_loci, drop=FALSE]
  if (ncol(phylo_gt)<2L) stop("Fewer than two polymorphic loci remain.", call.=FALSE)
  
  
  #### 2.3 Metadata groups and palettes
  
  clean_phylo_group <- function(x) { x <- as.character(x); x[is.na(x) | trimws(x)==""] <- "Missing"; x }
  site_by_sample <- setNames(clean_phylo_group(phylo_meta[, site_col_name, drop=TRUE]), phylo_meta$.phylo_sample)
  species_by_sample <- setNames(clean_phylo_group(phylo_meta[, species_col_name, drop=TRUE]), phylo_meta$.phylo_sample)
  
  complete_phylo_colours <- function(values, supplied_palette) {
    levels_needed <- unique(unname(values))
    if (is.null(supplied_palette) || is.null(names(supplied_palette))) supplied_palette <- character(0)
    supplied_palette <- supplied_palette[!is.na(names(supplied_palette)) & names(supplied_palette)!="" & !duplicated(names(supplied_palette))]
    missing_levels <- setdiff(levels_needed, names(supplied_palette))
    if (length(missing_levels)>0L) {
      new_colours <- colorspace::rainbow_hcl(length(missing_levels)); names(new_colours) <- missing_levels
      if ("Missing" %in% missing_levels) new_colours["Missing"] <- "grey50"
      supplied_palette <- c(supplied_palette, new_colours)
    }
    supplied_palette[levels_needed]
  }
  
  complete_phylo_shapes <- function(values, supplied_shapes) {
    levels_needed <- unique(unname(values))
    if (is.null(supplied_shapes) || is.null(names(supplied_shapes))) supplied_shapes <- integer(0)
    supplied_shapes <- supplied_shapes[!is.na(names(supplied_shapes)) & names(supplied_shapes)!="" & !duplicated(names(supplied_shapes))]
    missing_levels <- setdiff(levels_needed, names(supplied_shapes))
    if (length(missing_levels)>0L) {
      new_shapes <- rep(c(21:25,3,4,8,12), length.out=length(missing_levels)); names(new_shapes) <- missing_levels
      if ("Missing" %in% missing_levels) new_shapes["Missing"] <- 4
      supplied_shapes <- c(supplied_shapes, new_shapes)
    }
    supplied_shapes[levels_needed]
  }
  
  site_phylo_colours <- complete_phylo_colours(site_by_sample, site_colours)
  species_phylo_colours <- complete_phylo_colours(species_by_sample, sp_colours)
  site_phylo_shapes <- complete_phylo_shapes(site_by_sample, site_shapes)
  species_phylo_shapes <- complete_phylo_shapes(species_by_sample, sp_shapes)
  
  
  #### 2.4 Plot labels, dimensions, legends and graphics devices
  
  get_phylo_sizes <- function(n_tips, tanglegram=FALSE) {
    if (tanglegram) return(list(label=if (n_tips<=50) 0.90 else if (n_tips<=100) 0.72 else if (n_tips<=200) 0.60 else 0.52,
                                point=if (n_tips<=50) 1.50 else if (n_tips<=100) 1.25 else if (n_tips<=200) 1.00 else 0.90,
                                support=if (n_tips<=50) 0.78 else if (n_tips<=100) 0.64 else if (n_tips<=200) 0.52 else 0.44))
    list(label=if (n_tips<=50) 0.96 else if (n_tips<=100) 0.78 else if (n_tips<=200) 0.68 else 0.62,
         point=if (n_tips<=50) 1.55 else if (n_tips<=100) 1.30 else if (n_tips<=200) 1.10 else 0.95,
         support=if (n_tips<=50) 0.84 else if (n_tips<=100) 0.68 else if (n_tips<=200) 0.56 else 0.48)
  }
  
  make_phylo_tip_labels <- function(sample_names, label_mode=c("species","site","species_site")) {
    label_mode <- match.arg(label_mode); sample_names <- as.character(sample_names)
    species_labels <- clean_phylo_group(unname(species_by_sample[sample_names]))
    site_labels <- clean_phylo_group(unname(site_by_sample[sample_names]))
    switch(label_mode,
           species=paste(species_labels, sample_names, sep=phylo_tree_tip_label_separator),
           site=paste(site_labels, sample_names, sep=phylo_tree_tip_label_separator),
           species_site=paste(species_labels, site_labels, sample_names, sep=phylo_tree_tip_label_separator))
  }
  
  phylo_tip_label_description <- function(label_mode) switch(label_mode, species="species and sample accession", site="site and sample accession", species_site="species, site and sample accession")
  
  get_phylo_legend_layout <- function(group_labels) {
    group_labels <- unique(as.character(group_labels)); group_labels <- group_labels[!is.na(group_labels) & trimws(group_labels)!=""]; n_groups <- length(group_labels)
    if (n_groups==0L) return(list(columns=1L, rows=1L, fraction=0.07, cex=0.80))
    legend_cex <- if (n_groups>120) 0.68 else if (n_groups>80) 0.72 else if (n_groups>50) 0.76 else if (n_groups>30) 0.80 else 0.86
    preferred_columns <- if (n_groups>100) 6L else if (n_groups>50) 5L else if (n_groups>20) 4L else if (n_groups>8) 3L else if (n_groups>4) 2L else 1L
    device_width_inches <- unname(grDevices::dev.size("in")[1]); if (!is.finite(device_width_inches) || device_width_inches<=0) device_width_inches <- 30
    maximum_label_width <- max(graphics::strwidth(group_labels, units="inches", cex=legend_cex), na.rm=TRUE)
    if (!is.finite(maximum_label_width) || maximum_label_width<=0) maximum_label_width <- 1
    maximum_columns_by_width <- max(1L, floor(device_width_inches*0.95/(maximum_label_width+0.55)))
    legend_columns <- max(1L, min(n_groups, preferred_columns, as.integer(phylo_tree_legend_max_columns), maximum_columns_by_width))
    legend_rows <- ceiling(n_groups/legend_columns); legend_fraction <- min(0.20, max(0.07, 0.040+0.0050*legend_rows))
    list(columns=legend_columns, rows=legend_rows, fraction=legend_fraction, cex=legend_cex)
  }
  
  save_pdf_png <- function(file_stem, width_cm, height_cm, draw_fun) {
    pdf_file <- file.path(phylo_dir, paste0(file_stem, ".pdf")); png_file <- file.path(phylo_dir, paste0(file_stem, ".png"))
    grDevices::pdf(pdf_file, width=width_cm/2.54, height=height_cm/2.54, useDingbats=FALSE, onefile=TRUE)
    tryCatch(draw_fun(), finally=grDevices::dev.off())
    png_res <- if (max(width_cm, height_cm)>phylo_tree_large_figure_threshold_cm) phylo_tree_large_png_res else phylo_tree_png_res
    if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_png(filename=png_file, width=width_cm, height=height_cm, units="cm", res=png_res, background="white")
    else if (capabilities("cairo")) grDevices::png(filename=png_file, width=width_cm, height=height_cm, units="cm", res=png_res, type="cairo-png", bg="white")
    else grDevices::png(filename=png_file, width=width_cm, height=height_cm, units="cm", res=png_res, bg="white")
    tryCatch(draw_fun(), finally=grDevices::dev.off())
  }
  
  
  #### 2.5 Rooting, tree and bootstrap-support helpers
  
  root_phylo_for_display <- function(tree) {
    if (ape::is.rooted(tree)) tree <- ape::unroot(tree)
    outgroup_tips <- character(0)
    if (length(configured_outgroup_taxa)>0L) {
      absent_taxa <- setdiff(configured_outgroup_taxa, unique(unname(species_by_sample[tree$tip.label])))
      if (length(absent_taxa)>0L) warning("Configured outgroup taxon/taxa absent from this tree: ", paste(absent_taxa, collapse=", "), call.=FALSE)
      outgroup_tips <- intersect(tree$tip.label, names(species_by_sample)[species_by_sample %in% configured_outgroup_taxa])
    }
    if (length(outgroup_tips)>0L) {
      if (length(outgroup_tips)>=ape::Ntip(tree)) stop("All tree tips were selected as outgroup samples.", call.=FALSE)
      if (length(outgroup_tips)>1L && !isTRUE(tryCatch(ape::is.monophyletic(tree, outgroup_tips), error=function(e) FALSE))) stop("The specified outgroup samples are not monophyletic in this tree.", call.=FALSE)
      tree <- ape::root(tree, outgroup=outgroup_tips, resolve.root=TRUE); method <- paste0("Outgroup rooted using ", paste(configured_outgroup_taxa, collapse=", "))
    } else {
      if (length(configured_outgroup_taxa)>0L) warning("No configured outgroup samples were available; midpoint rooting was used instead.", call.=FALSE)
      tree <- phangorn::midpoint(tree); method <- "Midpoint rooted for visualisation only"
    }
    list(tree=tree, method=method, outgroup_tips=outgroup_tips)
  }
  
  make_phylo_nj_tree <- function(x) {
    x <- as.matrix(x); storage.mode(x) <- "numeric"; dm <- as.matrix(stats::dist(x, method="euclidean"))
    if (any(!is.finite(dm))) stop("A distance matrix contains NA or infinite values.", call.=FALSE)
    if (max(dm, na.rm=TRUE)<=0) stop("A distance matrix contains no positive genetic distances.", call.=FALSE)
    ape::unroot(ape::nj(stats::as.dist(dm)))
  }
  
  as_multi_phylo <- function(trees) { if (inherits(trees, "multiPhylo")) return(trees); trees <- as.list(trees); class(trees) <- "multiPhylo"; trees }
  
  prune_bootstrap_trees <- function(trees, tips) {
    result <- lapply(as_multi_phylo(trees), function(tree) ape::keep.tip(tree, intersect(tips, tree$tip.label))); class(result) <- "multiPhylo"; result
  }
  
  calculate_bootstrap_support <- function(tree, bootstrap_trees) {
    if (is.null(bootstrap_trees) || length(bootstrap_trees)<1L) return(rep(NA_real_, tree$Nnode))
    bootstrap_trees <- as_multi_phylo(bootstrap_trees); support <- 100*ape::prop.clades(tree, bootstrap_trees, rooted=FALSE)/length(bootstrap_trees)
    if (ape::is.rooted(tree)) {
      root_node <- setdiff(unique(tree$edge[,1]), tree$edge[,2]); root_index <- root_node-ape::Ntip(tree)
      if (length(root_index)==1L && root_index>=1L && root_index<=length(support)) support[root_index] <- NA_real_
    }
    support
  }
  
  attach_phylo_bootstrap_support <- function(tree, bootstrap_trees) {
    support <- calculate_bootstrap_support(tree, bootstrap_trees); tree$node.label <- ifelse(is.na(support), NA_character_, formatC(support, format="f", digits=1))
    attr(tree, "bootstrap_support") <- support; tree
  }
  
  extract_support_table <- function(tree, method) {
    support <- suppressWarnings(as.numeric(tree$node.label))
    data.frame(Method=method, Node=ape::Ntip(tree)+seq_along(support), Bootstrap_support=support, stringsAsFactors=FALSE)
  }
  
  prepare_supported_display_tree <- function(tree) {
    tree$node.label <- NULL; rooted <- root_phylo_for_display(tree); rooted$tree <- ape::ladderize(rooted$tree, right=FALSE); rooted
  }
  
  
  #### 2.6 Timing and checkpoint helpers
  
  format_phylo_duration <- function(seconds) {
    seconds <- max(0, round(as.numeric(seconds))); days <- seconds%/%86400; seconds <- seconds%%86400; hours <- seconds%/%3600; seconds <- seconds%%3600; minutes <- seconds%/%60; seconds <- seconds%%60
    if (days>0) sprintf("%dd %02dh %02dm %02ds", days, hours, minutes, seconds) else if (hours>0) sprintf("%02dh %02dm %02ds", hours, minutes, seconds) else sprintf("%02dm %02ds", minutes, seconds)
  }
  
  phylo_progress_message <- function(stage, i, total, durations, start_time, recent_n=5L) {
    completed_times <- durations[seq_len(i)]; completed_times <- completed_times[is.finite(completed_times)]
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs")); last_time <- if (length(completed_times)>0L) tail(completed_times,1) else NA_real_
    mean_time <- if (length(completed_times)>0L) mean(completed_times) else NA_real_; recent_time <- if (length(completed_times)>0L) mean(tail(completed_times, min(recent_n,length(completed_times)))) else NA_real_
    eta <- if (is.finite(recent_time)) recent_time*(total-i) else NA_real_; finish_time <- if (is.finite(eta)) Sys.time()+eta else as.POSIXct(NA)
    message(sprintf("[%s] %d/%d completed (%5.1f%%) | last: %s | mean: %s | elapsed: %s | ETA: %s | finish: %s", stage, i, total, 100*i/total,
                    if (is.finite(last_time)) format_phylo_duration(last_time) else "unknown", if (is.finite(mean_time)) format_phylo_duration(mean_time) else "unknown",
                    format_phylo_duration(elapsed), if (is.finite(eta)) format_phylo_duration(eta) else "unknown",
                    if (!is.na(finish_time)) format(finish_time, "%Y-%m-%d %H:%M:%S") else "unknown"))
    flush.console()
  }
  
  load_phylo_checkpoint <- function(checkpoint_file, total_reps, seed, resume=TRUE) {
    if (resume && !is.null(checkpoint_file) && file.exists(checkpoint_file)) {
      checkpoint <- readRDS(checkpoint_file); trees <- checkpoint$trees; durations <- checkpoint$durations
      if (length(trees)>total_reps) trees <- trees[seq_len(total_reps)]
      if (length(durations)>total_reps) durations <- durations[seq_len(total_reps)]
      if (!is.null(checkpoint$rng_state)) assign(".Random.seed", checkpoint$rng_state, envir=.GlobalEnv)
      message("[Checkpoint] Resuming from ", length(trees), " completed replicate(s): ", checkpoint_file)
    } else { set.seed(seed); trees <- list(); durations <- numeric(0) }
    list(trees=trees, durations=durations)
  }
  
  save_phylo_checkpoint <- function(trees, durations, checkpoint_file, total_reps) {
    if (is.null(checkpoint_file)) return(invisible(NULL))
    saveRDS(list(trees=trees, durations=durations, target_replicates=total_reps, rng_state=get0(".Random.seed", envir=.GlobalEnv, inherits=FALSE), saved_at=Sys.time()), checkpoint_file)
    invisible(NULL)
  }
  
  
  #### 2.7 Parallel and bootstrap functions

  make_phylo_rng_streams <- function(seed, n) {
    old_kind <- RNGkind(); on.exit(do.call(RNGkind, as.list(old_kind)), add=TRUE)
    RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
    streams <- vector("list",n); current <- .Random.seed
    for (i in seq_len(n)) { streams[[i]] <- current; current <- parallel::nextRNGStream(current) }
    streams
  }

  phylo_lapply <- function(X, FUN, use_parallel=TRUE, ...) {
    worker_count <- min(length(X),phylo_n_cores)
    if (isTRUE(phylo_use_parallel) && isTRUE(use_parallel) && worker_count>1L && .Platform$OS.type!="windows") {
      parallel::mclapply(X,FUN,...,mc.cores=worker_count,mc.preschedule=FALSE,mc.set.seed=TRUE)
    } else lapply(X,FUN,...)
  }

  phylo_parallel_progress <- function(stage, completed, total, completed_at_start, run_start, active_workers) {
    elapsed <- as.numeric(difftime(Sys.time(),run_start,units="secs")); completed_now <- completed-completed_at_start
    throughput <- if (elapsed>0 && completed_now>0) completed_now/elapsed else NA_real_
    eta <- if (is.finite(throughput) && throughput>0) (total-completed)/throughput else NA_real_
    finish_time <- if (is.finite(eta)) Sys.time()+eta else as.POSIXct(NA)
    message(sprintf("[%s] %d/%d completed (%5.1f%%) | workers: %d | current-run elapsed: %s | ETA: %s | finish: %s", stage,completed,total,100*completed/total,
                    active_workers,format_phylo_duration(elapsed),if (is.finite(eta)) format_phylo_duration(eta) else "unknown",
                    if (!is.na(finish_time)) format(finish_time,"%Y-%m-%d %H:%M:%S") else "unknown")); flush.console()
  }

  run_nj_bootstrap_progress <- function(x, B, FUN, seed, checkpoint_file=NULL, progress_every=1L, checkpoint_every=5L, resume=TRUE, recent_n=5L) {
    state <- load_phylo_checkpoint(checkpoint_file,B,seed,resume); trees <- state$trees; durations <- state$durations; first_rep <- length(trees)+1L
    completed_at_start <- length(trees); run_start <- Sys.time(); rng_streams <- make_phylo_rng_streams(seed,B)
    message(sprintf("[NJ bootstrap] Starting %d locus-bootstrap replicates using %d samples, %s loci and %d worker core(s).",B,nrow(x),format(ncol(x),big.mark=","),phylo_n_cores))
    if (first_rep>B) message("[NJ bootstrap] All requested replicates were already present in the checkpoint.")

    if (first_rep<=B) {
      pending <- seq.int(first_rep,B); batch_size <- if (isTRUE(phylo_parallel_bootstraps)) max(1L,phylo_n_cores) else 1L
      batches <- split(pending,ceiling(seq_along(pending)/batch_size))
      for (batch_ids in batches) {
        results <- phylo_lapply(batch_ids,function(i) {
          replicate_start <- Sys.time(); assign(".Random.seed",rng_streams[[i]],envir=.GlobalEnv)
          tryCatch({
            sampled_loci <- sample.int(ncol(x),ncol(x),replace=TRUE)
            list(ok=TRUE,tree=FUN(x[,sampled_loci,drop=FALSE]),duration=as.numeric(difftime(Sys.time(),replicate_start,units="secs")),error=NA_character_)
          },error=function(e) list(ok=FALSE,tree=NULL,duration=as.numeric(difftime(Sys.time(),replicate_start,units="secs")),error=conditionMessage(e)))
        },use_parallel=phylo_parallel_bootstraps)

        result_ok <- vapply(results,function(z) is.list(z) && isTRUE(z$ok),logical(1))
        if (any(!result_ok)) stop("NJ bootstrap replicate(s) failed: ",paste(vapply(results[!result_ok],function(z) if (is.list(z)) z$error else as.character(z),character(1)),collapse=" | "),call.=FALSE)
        for (j in seq_along(batch_ids)) { i <- batch_ids[j]; trees[[i]] <- results[[j]]$tree; durations[i] <- results[[j]]$duration }

        completed <- max(batch_ids); phylo_parallel_progress("NJ bootstrap",completed,B,completed_at_start,run_start,min(length(batch_ids),phylo_n_cores))
        if (!is.null(checkpoint_file)) { save_phylo_checkpoint(trees,durations,checkpoint_file,B); message("[NJ bootstrap] Checkpoint saved after replicate ",completed,".") }
      }
    }

    trees <- trees[seq_len(B)]; class(trees) <- "multiPhylo"
    if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) file.remove(checkpoint_file)
    message("[NJ bootstrap] Completed ",B," replicates; current server run wall time: ",format_phylo_duration(difftime(Sys.time(),run_start,units="secs")),"."); trees
  }

  run_ml_bootstrap_progress <- function(fit, B, seed, checkpoint_file=NULL, progress_every=1L, checkpoint_every=5L, resume=TRUE, recent_n=5L) {
    state <- load_phylo_checkpoint(checkpoint_file,B,seed,resume); trees <- state$trees; durations <- state$durations; first_rep <- length(trees)+1L
    completed_at_start <- length(trees); run_start <- Sys.time(); rng_streams <- make_phylo_rng_streams(seed,B)
    message(sprintf("[ML bootstrap] Starting %d replicates under %s+ASC using %d worker core(s).",B,phylo_ml_model,phylo_n_cores))
    message("[ML bootstrap] Each worker fits one independent bootstrap replicate; checkpoints are saved after every completed batch.")
    if (first_rep>B) message("[ML bootstrap] All requested replicates were already present in the checkpoint.")

    if (first_rep<=B) {
      pending <- seq.int(first_rep,B); batch_size <- if (isTRUE(phylo_parallel_bootstraps)) max(1L,phylo_n_cores) else 1L
      batches <- split(pending,ceiling(seq_along(pending)/batch_size))
      for (batch_ids in batches) {
        message(sprintf("[ML bootstrap] Starting replicate batch %d-%d at %s.",min(batch_ids),max(batch_ids),format(Sys.time(),"%Y-%m-%d %H:%M:%S"))); flush.console()
        results <- phylo_lapply(batch_ids,function(i) {
          replicate_start <- Sys.time(); assign(".Random.seed",rng_streams[[i]],envir=.GlobalEnv)
          tryCatch({
            one_bootstrap <- phangorn::bootstrap.pml(fit,bs=1L,trees=TRUE,multicore=FALSE,optNni=TRUE,control=phangorn::pml.control(trace=0))
            list(ok=TRUE,tree=one_bootstrap[[1]],duration=as.numeric(difftime(Sys.time(),replicate_start,units="secs")),error=NA_character_)
          },error=function(e) list(ok=FALSE,tree=NULL,duration=as.numeric(difftime(Sys.time(),replicate_start,units="secs")),error=conditionMessage(e)))
        },use_parallel=phylo_parallel_bootstraps)

        result_ok <- vapply(results,function(z) is.list(z) && isTRUE(z$ok),logical(1))
        if (any(!result_ok)) stop("ML bootstrap replicate(s) failed: ",paste(vapply(results[!result_ok],function(z) if (is.list(z)) z$error else as.character(z),character(1)),collapse=" | "),call.=FALSE)
        for (j in seq_along(batch_ids)) { i <- batch_ids[j]; trees[[i]] <- results[[j]]$tree; durations[i] <- results[[j]]$duration }

        completed <- max(batch_ids); phylo_parallel_progress("ML bootstrap",completed,B,completed_at_start,run_start,min(length(batch_ids),phylo_n_cores))
        if (!is.null(checkpoint_file)) { save_phylo_checkpoint(trees,durations,checkpoint_file,B); message("[ML bootstrap] Checkpoint saved after replicate ",completed,".") }
      }
    }

    trees <- trees[seq_len(B)]; class(trees) <- "multiPhylo"
    if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) file.remove(checkpoint_file)
    message("[ML bootstrap] Completed ",B," replicates; current server run wall time: ",format_phylo_duration(difftime(Sys.time(),run_start,units="secs")),"."); trees
  }


  #### 2.8 Stand-alone NJ/ML plotting
  
  draw_supported_phylo_tree <- function(tree, bootstrap_trees, group_by_sample, colour_palette, shape_palette, group_title, label_mode, method_title) {
    label_mode <- match.arg(label_mode, c("species","site","species_site")); display_rooting <- prepare_supported_display_tree(tree); display_tree <- display_rooting$tree
    support <- calculate_bootstrap_support(display_tree, bootstrap_trees); sample_names <- as.character(display_tree$tip.label); tip_groups <- clean_phylo_group(unname(group_by_sample[sample_names]))
    plot_tree <- display_tree; plot_tree$tip.label <- make_phylo_tip_labels(sample_names, label_mode)
    tip_colours <- unname(colour_palette[tip_groups]); tip_shapes <- unname(shape_palette[tip_groups]); tip_colours[is.na(tip_colours)] <- "grey50"; tip_shapes[is.na(tip_shapes)] <- 4
    tip_outlines <- ifelse(tip_shapes %in% 21:25, "black", tip_colours); label_colours <- if (phylo_tree_colour_sample_labels) tip_colours else rep("black", length(tip_colours))
    sizes <- get_phylo_sizes(ape::Ntip(plot_tree))
    tree_span <- if (!phylo_tree_show_branch_lengths || is.null(plot_tree$edge.length)) 1 else max(ape::node.depth.edgelength(plot_tree), na.rm=TRUE)
    if (!is.finite(tree_span) || tree_span<=0) tree_span <- 1
    present_groups <- unique(tip_groups); legend_colours <- unname(colour_palette[present_groups]); legend_shapes <- unname(shape_palette[present_groups])
    legend_colours[is.na(legend_colours)] <- "grey50"; legend_shapes[is.na(legend_shapes)] <- 4; legend_outlines <- ifelse(legend_shapes %in% 21:25, "black", legend_colours)
    legend_layout <- get_phylo_legend_layout(present_groups)
    graphics::layout(matrix(c(1,2), ncol=1), heights=c(1-legend_layout$fraction, legend_layout$fraction))
    graphics::par(mar=c(2.2,2,5.0,3), xpd=NA)
    ape::plot.phylo(plot_tree, type="phylogram", direction="rightwards", use.edge.length=phylo_tree_show_branch_lengths, show.tip.label=TRUE, tip.color=label_colours,
                    cex=sizes$label, edge.width=0.70, label.offset=tree_span*0.008, align.tip.label=FALSE, no.margin=FALSE)
    ape::tiplabels(pch=tip_shapes, col=tip_outlines, bg=tip_colours, cex=sizes$point, lwd=0.7)
    supported_nodes <- which(!is.na(support) & support>=phylo_support_to_plot)
    if (length(supported_nodes)>0L) ape::nodelabels(text=format(round(support[supported_nodes]), trim=TRUE), node=ape::Ntip(plot_tree)+supported_nodes, frame="none", cex=sizes$support, font=2, col="black")
    graphics::title(main=paste0(species, ": ", method_title, "\nTip labels show ", phylo_tip_label_description(label_mode),
                                "\nNode labels are locus-bootstrap support (%) | ", display_rooting$method), cex.main=1.25)
    graphics::par(mar=c(0.2,2,0.2,2), xpd=NA); graphics::plot.new()
    graphics::legend("top", inset=c(0,0.01), legend=present_groups, title=group_title, pch=legend_shapes, col=legend_outlines, pt.bg=legend_colours,
                     pt.cex=1.20, cex=legend_layout$cex, bty="n", ncol=legend_layout$columns, x.intersp=0.45, y.intersp=1)
    invisible(list(display_tree=display_tree, plot_tree=plot_tree, support=support, rooting_method=display_rooting$method))
  }
  
  save_supported_phylo_tree <- function(tree, bootstrap_trees, group_by_sample, colour_palette, shape_palette, group_title, label_mode, method_title, file_stem) {
    label_mode <- match.arg(label_mode, c("species","site","species_site")); sample_names <- as.character(tree$tip.label)
    combined_tip_labels <- make_phylo_tip_labels(sample_names, label_mode); present_groups <- unique(clean_phylo_group(unname(group_by_sample[sample_names])))
    n_tips <- ape::Ntip(tree); max_tip_chars <- max(nchar(combined_tip_labels), na.rm=TRUE); max_legend_chars <- max(nchar(present_groups), na.rm=TRUE)
    width_cm <- min(phylo_tree_max_width_cm, max(phylo_tree_min_width_cm, 45+max_tip_chars*1.00, 50+max_legend_chars*0.70))
    height_cm <- min(phylo_tree_max_height_cm, max(phylo_tree_min_height_cm, 20+n_tips*phylo_tree_height_per_tip_cm))
    message("[Tree output] Saving ", file_stem, " at ", round(width_cm,1), " × ", round(height_cm,1), " cm.")
    save_pdf_png(file_stem, width_cm, height_cm, function() draw_supported_phylo_tree(tree, bootstrap_trees, group_by_sample, colour_palette, shape_palette, group_title, label_mode, method_title))
  }
  
  
  #### 2.8B COLLAPSED MONOPHYLETIC-SPECIES TREES #########################################
  
  # Identify species whose samples form one monophyletic clade. Monophyletic species with
  # more than one sample are reduced to one representative plotting tip. Non-monophyletic
  # species remain expanded so their individual sample placements can still be examined.
  collapse_monophyletic_species_tree <- function(tree, species_lookup) {
    sample_names <- as.character(tree$tip.label)
    species_labels <- clean_phylo_group(unname(species_lookup[sample_names]))
    species_split <- split(sample_names, species_labels)
    
    replacement_labels <- setNames(sample_names, sample_names)
    kept_groups <- setNames(species_labels, sample_names)
    collapsed_by_tip <- setNames(rep(FALSE, length(sample_names)), sample_names)
    collapsed_n_by_tip <- setNames(rep(1L, length(sample_names)), sample_names)
    tips_to_drop <- character(0)
    
    collapse_table <- data.frame(
      Species=names(species_split), N_samples=integer(length(species_split)), Monophyletic=logical(length(species_split)),
      Collapsed_in_plot=logical(length(species_split)), Representative_tip=NA_character_, Plot_label=NA_character_,
      stringsAsFactors=FALSE
    )
    
    for (i in seq_along(species_split)) {
      sp <- names(species_split)[i]
      tips <- species_split[[i]]
      is_mono <- length(tips)==1L || isTRUE(tryCatch(ape::is.monophyletic(tree, tips), error=function(e) FALSE))
      
      collapse_table$N_samples[i] <- length(tips)
      collapse_table$Monophyletic[i] <- is_mono
      
      if (length(tips)==1L) {
        replacement_labels[tips] <- sp
        collapse_table$Collapsed_in_plot[i] <- FALSE
        collapse_table$Representative_tip[i] <- tips
        collapse_table$Plot_label[i] <- sp
        
      } else if (is_mono) {
        representative_tip <- tips[1]
        collapsed_label <- paste0(sp, " (n=", length(tips), ")")
        
        replacement_labels[representative_tip] <- collapsed_label
        collapsed_by_tip[representative_tip] <- TRUE
        collapsed_n_by_tip[representative_tip] <- length(tips)
        tips_to_drop <- c(tips_to_drop, tips[-1])
        
        collapse_table$Collapsed_in_plot[i] <- TRUE
        collapse_table$Representative_tip[i] <- representative_tip
        collapse_table$Plot_label[i] <- collapsed_label
        
      } else {
        replacement_labels[tips] <- paste(sp, tips, sep=phylo_tree_tip_label_separator)
        collapse_table$Collapsed_in_plot[i] <- FALSE
        collapse_table$Representative_tip[i] <- NA_character_
        collapse_table$Plot_label[i] <- "Individual samples retained because the species is not monophyletic"
      }
    }
    
    collapsed_tree <- if (length(tips_to_drop)>0L) ape::drop.tip(tree, tips_to_drop, trim.internal=TRUE, subtree=FALSE) else tree
    retained_original_tips <- as.character(collapsed_tree$tip.label)
    
    list(
      tree=collapsed_tree,
      original_tip_names=retained_original_tips,
      plot_labels=unname(replacement_labels[retained_original_tips]),
      tip_groups=unname(kept_groups[retained_original_tips]),
      collapsed_tip=unname(collapsed_by_tip[retained_original_tips]),
      collapsed_n=unname(collapsed_n_by_tip[retained_original_tips]),
      collapse_table=collapse_table[order(collapse_table$Species), , drop=FALSE]
    )
  }
  
  
  # Root the display copy, collapse monophyletic species, prune bootstrap trees to the
  # retained representative tips and recalculate support for the reduced plotting tree.
  prepare_collapsed_species_plot_data <- function(tree, bootstrap_trees) {
    display_rooting <- prepare_supported_display_tree(tree)
    collapsed <- collapse_monophyletic_species_tree(display_rooting$tree, species_by_sample)
    
    bootstrap_common <- prune_bootstrap_trees(bootstrap_trees, collapsed$original_tip_names)
    collapsed_support <- calculate_bootstrap_support(collapsed$tree, bootstrap_common)
    
    plot_tree <- collapsed$tree
    plot_tree$tip.label <- collapsed$plot_labels
    
    list(
      plot_tree=plot_tree,
      support=collapsed_support,
      tip_groups=collapsed$tip_groups,
      collapsed_tip=collapsed$collapsed_tip,
      collapsed_n=collapsed$collapsed_n,
      collapse_table=collapsed$collapse_table,
      rooting_method=display_rooting$method
    )
  }
  
  
  # Draw collapsed monophyletic species as large sideways triangles. The narrow point
  # begins at the end of the retained branch and the triangle widens to the right.
  draw_collapsed_species_phylo_tree <- function(prepared, method_title) {
    plot_tree <- prepared$plot_tree
    support <- prepared$support
    tip_groups <- clean_phylo_group(prepared$tip_groups)
    collapsed_tip <- as.logical(prepared$collapsed_tip)
    
    tip_colours <- unname(species_phylo_colours[tip_groups])
    tip_shapes <- unname(species_phylo_shapes[tip_groups])
    tip_colours[is.na(tip_colours)] <- "grey50"
    tip_shapes[is.na(tip_shapes)] <- 4
    
    tip_outlines <- ifelse(tip_shapes %in% 21:25, "black", tip_colours)
    label_colours <- if (phylo_tree_colour_sample_labels) tip_colours else rep("black", length(tip_colours))
    
    sizes <- get_phylo_sizes(ape::Ntip(plot_tree))
    label_cex <- sizes$label * phylo_collapsed_label_scale
    
    tree_span <- if (phylo_tree_show_branch_lengths && !is.null(plot_tree$edge.length)) {
      max(ape::node.depth.edgelength(plot_tree), na.rm=TRUE)
    } else {
      max(ape::node.depth(plot_tree), na.rm=TRUE)
    }
    if (!is.finite(tree_span) || tree_span<=0) tree_span <- 1
    
    present_groups <- unique(tip_groups)
    legend_colours <- unname(species_phylo_colours[present_groups])
    legend_shapes <- unname(species_phylo_shapes[present_groups])
    legend_colours[is.na(legend_colours)] <- "grey50"
    legend_shapes[is.na(legend_shapes)] <- 4
    legend_outlines <- ifelse(legend_shapes %in% 21:25, "black", legend_colours)
    legend_layout <- get_phylo_legend_layout(present_groups)
    
    graphics::layout(matrix(c(1,2), ncol=1), heights=c(1-legend_layout$fraction, legend_layout$fraction))
    graphics::par(mar=c(2.2,2,5.4,6), xpd=NA)
    
    # Invisible labels make plot.phylo reserve sufficient space. Labels are redrawn
    # manually after the triangles so they begin beyond the broad end of each triangle.
    invisible_label_colour <- rep(grDevices::adjustcolor("black", alpha.f=0), ape::Ntip(plot_tree))
    reserve_offset <- tree_span * (phylo_collapsed_triangle_length_fraction + 0.018)
    
    ape::plot.phylo(
      plot_tree, type="phylogram", direction="rightwards", use.edge.length=phylo_tree_show_branch_lengths,
      show.tip.label=TRUE, tip.color=invisible_label_colour, cex=label_cex, edge.width=0.70,
      label.offset=reserve_offset, align.tip.label=FALSE, no.margin=FALSE
    )
    
    ape_plot_environment <- getFromNamespace(".PlotPhyloEnv", "ape")
    last_plot <- get0("last_plot.phylo", envir=ape_plot_environment, inherits=FALSE)
    if (is.null(last_plot) || is.null(last_plot$xx) || is.null(last_plot$yy)) stop("Tip coordinates could not be recovered after plotting the collapsed phylogeny.", call.=FALSE)
    
    tip_x <- last_plot$xx[seq_len(ape::Ntip(plot_tree))]
    tip_y <- last_plot$yy[seq_len(ape::Ntip(plot_tree))]
    plot_x_span <- diff(graphics::par("usr")[1:2])
    
    
    triangle_length <- plot_x_span * phylo_collapsed_triangle_length_fraction
    triangle_half_height <- phylo_collapsed_triangle_half_height
    label_gap <- plot_x_span * 0.008
    
    # Ordinary symbols are retained for singleton species and individual samples from
    # species that are not monophyletic.
    ordinary_tips <- which(!collapsed_tip)
    if (length(ordinary_tips)>0L) {
      ape::tiplabels(
        tip=ordinary_tips, pch=tip_shapes[ordinary_tips], col=tip_outlines[ordinary_tips],
        bg=tip_colours[ordinary_tips], cex=sizes$point, lwd=0.7
      )
    }
    
    # Draw one large, right-expanding triangle for every collapsed species clade.
    collapsed_tips <- which(collapsed_tip)
    if (length(collapsed_tips)>0L) {
      for (i in collapsed_tips) {
        graphics::polygon(
          x=c(tip_x[i], tip_x[i]+triangle_length, tip_x[i]+triangle_length),
          y=c(tip_y[i], tip_y[i]-triangle_half_height, tip_y[i]+triangle_half_height),
          col=tip_colours[i], border="black", lwd=phylo_collapsed_triangle_border_width
        )
      }
    }
    
    # Draw ordinary and collapsed labels separately so collapsed labels are positioned
    # after the broad end of the triangle.
    if (length(ordinary_tips)>0L) {
      graphics::text(
        x=tip_x[ordinary_tips]+label_gap, y=tip_y[ordinary_tips], labels=plot_tree$tip.label[ordinary_tips],
        pos=4, offset=0, cex=label_cex, col=label_colours[ordinary_tips]
      )
    }
    
    if (length(collapsed_tips)>0L) {
      graphics::text(
        x=tip_x[collapsed_tips]+triangle_length+label_gap, y=tip_y[collapsed_tips], labels=plot_tree$tip.label[collapsed_tips],
        pos=4, offset=0, cex=label_cex, col=label_colours[collapsed_tips], font=2
      )
    }
    
    supported_nodes <- which(!is.na(support) & support>=phylo_support_to_plot)
    if (length(supported_nodes)>0L) {
      ape::nodelabels(
        text=format(round(support[supported_nodes]), trim=TRUE), node=ape::Ntip(plot_tree)+supported_nodes,
        frame="none", cex=sizes$support, font=2, col="black"
      )
    }
    
    graphics::title(
      main=paste0(
        species, ": ", method_title,
        "\nLarge coloured triangles represent monophyletic species collapsed to one branch; n is the number of samples represented.",
        "\nNon-monophyletic species remain expanded | Node labels are locus-bootstrap support (%) | ", prepared$rooting_method
      ),
      cex.main=1.20
    )
    
    graphics::par(mar=c(0.2,2,0.2,2), xpd=NA)
    graphics::plot.new()
    graphics::legend(
      "top", inset=c(0,0.01), legend=present_groups, title=species_col_name, pch=legend_shapes,
      col=legend_outlines, pt.bg=legend_colours, pt.cex=1.20, cex=legend_layout$cex,
      bty="n", ncol=legend_layout$columns, x.intersp=0.45, y.intersp=1
    )
    
    invisible(prepared)
  }
  
  
  save_collapsed_species_phylo_tree <- function(tree, bootstrap_trees, method_title, file_stem) {
    prepared <- prepare_collapsed_species_plot_data(tree, bootstrap_trees)
    
    plot_labels <- as.character(prepared$plot_tree$tip.label)
    present_groups <- unique(clean_phylo_group(prepared$tip_groups))
    n_tips <- ape::Ntip(prepared$plot_tree)
    
    max_tip_chars <- max(nchar(plot_labels), na.rm=TRUE)
    max_legend_chars <- max(nchar(present_groups), na.rm=TRUE)
    
    width_cm <- min(phylo_tree_max_width_cm, max(phylo_tree_min_width_cm, 55+max_tip_chars*1.00, 50+max_legend_chars*0.70))
    height_cm <- min(phylo_tree_max_height_cm, max(45, 20+n_tips*phylo_tree_height_per_tip_cm))
    
    message("[Collapsed species tree] Saving ", file_stem, " at ", round(width_cm,1), " × ", round(height_cm,1), " cm.")
    
    utils::write.csv(prepared$collapse_table, file.path(phylo_dir, paste0(file_stem, "_collapse_summary.csv")), row.names=FALSE)
    save_pdf_png(file_stem, width_cm, height_cm, function() draw_collapsed_species_phylo_tree(prepared, method_title))
  }
  
  #### 2.9 Result objects and replot-only loader
  
  nj_tree <- nj_supported_tree <- nj_bootstrap_trees <- NULL
  ml_tree <- ml_fit <- ml_supported_tree <- ml_bootstrap_trees <- NULL
  prepared_tanglegram <- NULL
  
  read_required_phylo_rds <- function(path, description) {
    if (!file.exists(path)) stop("Replot-only mode requires ", description, ", but the file was not found: ", path, call.=FALSE)
    readRDS(path)
  }
  
  if (phylo_replot_only) {
    message("[Replot-only] Loading previously saved tree and bootstrap objects from: ", phylo_dir)
    if (phylo_run_nj) {
      nj_supported_tree <- read_required_phylo_rds(file.path(phylo_dir, "Bootstrap_NJ_tree_supported.rds"), "the supported NJ tree")
      nj_bootstrap_trees <- read_required_phylo_rds(file.path(phylo_dir, "Bootstrap_NJ_replicate_trees.rds"), "the NJ bootstrap replicate trees")
    }
    if (phylo_run_ml) {
      ml_supported_tree <- read_required_phylo_rds(file.path(phylo_dir, paste0("Bootstrap_ML_", phylo_ml_model, "_ASC_tree_supported.rds")), "the supported ML tree")
      ml_bootstrap_trees <- read_required_phylo_rds(file.path(phylo_dir, "Bootstrap_ML_replicate_trees.rds"), "the ML bootstrap replicate trees")
      ml_fit_file <- file.path(phylo_dir, paste0("Bootstrap_ML_", phylo_ml_model, "_ASC_model_fit.rds"))
      if (file.exists(ml_fit_file)) ml_fit <- readRDS(ml_fit_file)
    }
  }
  
  
  #### 3. NEIGHBOUR-JOINING ANALYSIS #####################################################
  
  if (phylo_run_nj && !phylo_replot_only) {
    nj_stage_start <- Sys.time(); message(sprintf("[NJ tree] Constructing the reference neighbour-joining tree at %s.", format(nj_stage_start, "%Y-%m-%d %H:%M:%S")))
    nj_tree <- make_phylo_nj_tree(phylo_gt)
    message(sprintf("[NJ tree] Reference tree constructed in %s | tips: %d | internal nodes: %d.", format_phylo_duration(difftime(Sys.time(), nj_stage_start, units="secs")), ape::Ntip(nj_tree), nj_tree$Nnode))
    
    nj_checkpoint_file <- file.path(phylo_dir, "NJ_bootstrap_checkpoint.rds")
    if (isTRUE(phylo_reset_nj_checkpoint) && file.exists(nj_checkpoint_file)) { message("[NJ bootstrap] Removing the existing checkpoint."); file.remove(nj_checkpoint_file) }
    
    nj_bootstrap_trees <- run_nj_bootstrap_progress(phylo_gt, phylo_nj_bootstrap_reps, make_phylo_nj_tree, phylo_bootstrap_seed, nj_checkpoint_file,
                                                    phylo_progress_every, phylo_checkpoint_every, phylo_resume_bootstraps, phylo_progress_recent_n)
    message("[NJ support] Calculating locus-bootstrap support values."); support_start <- Sys.time()
    nj_supported_tree <- attach_phylo_bootstrap_support(nj_tree, nj_bootstrap_trees)
    message(sprintf("[NJ support] Support calculation completed in %s.", format_phylo_duration(difftime(Sys.time(), support_start, units="secs"))))
    
    ape::write.tree(nj_supported_tree, file=file.path(phylo_dir, "Bootstrap_NJ_tree_supported.newick"))
    saveRDS(nj_supported_tree, file=file.path(phylo_dir, "Bootstrap_NJ_tree_supported.rds"))
    if (phylo_save_bootstrap_trees) saveRDS(nj_bootstrap_trees, file=file.path(phylo_dir, "Bootstrap_NJ_replicate_trees.rds"))
    message(sprintf("[NJ complete] NJ analysis completed in %s.", format_phylo_duration(difftime(Sys.time(), nj_stage_start, units="secs"))))
  }
  
  
  #### 4. MAXIMUM-LIKELIHOOD ANALYSIS ####################################################
  
  if (phylo_run_ml && !phylo_replot_only) {
    ml_stage_start <- Sys.time(); message(sprintf("[ML setup] Preparing the genotype-state ML analysis at %s.", format(ml_stage_start, "%Y-%m-%d %H:%M:%S")))
    
    ml_character_matrix <- matrix(as.character(phylo_gt), nrow=nrow(phylo_gt), ncol=ncol(phylo_gt), dimnames=dimnames(phylo_gt))
    ml_character_matrix[is.na(ml_character_matrix)] <- "?"
    ml_contrast <- matrix(c(1,0,0, 0,1,0, 0,0,1, 1,1,1), nrow=4, byrow=TRUE, dimnames=list(c("0","1","2","?"), c("0","1","2")))
    ml_phyDat <- phangorn::phyDat(ml_character_matrix, type="USER", contrast=ml_contrast)
    ml_start_tree <- if (!is.null(nj_tree)) nj_tree else make_phylo_nj_tree(phylo_gt)
    ml_start_tree <- ape::unroot(ml_start_tree); ml_start_tree$node.label <- NULL
    message(sprintf("[ML setup] Three-state phyDat prepared: %d taxa × %s polymorphic loci.", nrow(phylo_gt), format(ncol(phylo_gt), big.mark=",")))
    
    optimise_fixed_ml_edges <- function(tree, start_length) {
      test_tree <- tree; test_tree$edge.length <- rep(start_length, nrow(test_tree$edge))
      test_fit <- phangorn::pml(test_tree, ml_phyDat, bf="equal", model=phylo_ml_model, ASC=TRUE)
      phangorn::optim.pml(test_fit, model=phylo_ml_model, optEdge=TRUE, optNni=FALSE, rearrangement="none",
                          control=phangorn::pml.control(epsilon=phylo_ml_edge_epsilon, maxit=phylo_ml_edge_maxit, trace=0))
    }
    
    message("[ML optimisation] Stage 1/3: optimising branch lengths on the fixed NJ starting topology.")
    ml_initial_edge_fits <- phylo_lapply(phylo_ml_edge_start_values, function(x) optimise_fixed_ml_edges(ml_start_tree,x), use_parallel=phylo_parallel_edge_starts)
    ml_initial_logLik <- vapply(ml_initial_edge_fits, function(x) as.numeric(x$logLik), numeric(1)); ml_initial_best_i <- which.max(ml_initial_logLik)
    ml_fit_start <- ml_initial_edge_fits[[ml_initial_best_i]]
    message(sprintf("[ML optimisation] Best initial fit: start edge %.4g | log-likelihood %.6f | maximum edge %.6f.", phylo_ml_edge_start_values[ml_initial_best_i], as.numeric(ml_fit_start$logLik), max(ml_fit_start$tree$edge.length)))
    
    message(sprintf("[ML optimisation] Stage 2/3: starting %s topology search. This optim.pml topology search uses one R process; parallel workers are used for the independent edge-start fits and bootstrap replicates.",phylo_ml_rearrangement)); ml_topology_start <- Sys.time()
    ml_fit_topology <- phangorn::optim.pml(ml_fit_start, model=phylo_ml_model, optEdge=TRUE, optNni=identical(phylo_ml_rearrangement,"NNI"),
                                           rearrangement=phylo_ml_rearrangement, control=phangorn::pml.control(trace=phylo_ml_optim_trace))
    message(sprintf("[ML optimisation] Topology search completed in %s | log-likelihood %.6f.", format_phylo_duration(difftime(Sys.time(), ml_topology_start, units="secs")), as.numeric(ml_fit_topology$logLik)))
    
    message("[ML optimisation] Stage 3/3: final multi-start fixed-topology branch-length refinement.")
    ml_final_edge_fits <- phylo_lapply(phylo_ml_edge_start_values, function(x) optimise_fixed_ml_edges(ml_fit_topology$tree,x), use_parallel=phylo_parallel_edge_starts)
    ml_final_logLik <- vapply(ml_final_edge_fits, function(x) as.numeric(x$logLik), numeric(1)); ml_final_best_i <- which.max(ml_final_logLik)
    ml_fit_edge_refined <- ml_final_edge_fits[[ml_final_best_i]]
    topology_logLik <- as.numeric(ml_fit_topology$logLik); refined_logLik <- as.numeric(ml_fit_edge_refined$logLik)
    
    if (refined_logLik>=topology_logLik-phylo_ml_loglik_tolerance) { ml_fit <- ml_fit_edge_refined; ml_selected_name <- "edge_refined" }
    else {
      ml_fit <- ml_fit_topology; ml_selected_name <- "topology_search"
      warning("The final fixed-topology edge refinement produced a lower likelihood than the topology-search fit. The higher-likelihood topology-search fit was retained; inspect ML_staged_optimisation_diagnostics.csv.", call.=FALSE)
    }
    
    message(sprintf("[ML optimisation] Final fit selected: %s | log-likelihood %.6f | maximum edge %.6f | total tree length %.6f.",
                    ml_selected_name, as.numeric(ml_fit$logLik), max(ml_fit$tree$edge.length), sum(ml_fit$tree$edge.length)))
    
    ml_optimisation_diagnostics <- data.frame(
      Stage=c(paste0("Initial fixed topology: ", phylo_ml_edge_start_values), "Topology search", paste0("Final fixed topology: ", phylo_ml_edge_start_values)),
      Starting_edge_length=c(phylo_ml_edge_start_values, NA_real_, phylo_ml_edge_start_values),
      Log_likelihood=c(ml_initial_logLik, topology_logLik, ml_final_logLik),
      Maximum_edge=c(vapply(ml_initial_edge_fits, function(x) max(x$tree$edge.length), numeric(1)), max(ml_fit_topology$tree$edge.length),
                     vapply(ml_final_edge_fits, function(x) max(x$tree$edge.length), numeric(1))),
      Total_tree_length=c(vapply(ml_initial_edge_fits, function(x) sum(x$tree$edge.length), numeric(1)), sum(ml_fit_topology$tree$edge.length),
                          vapply(ml_final_edge_fits, function(x) sum(x$tree$edge.length), numeric(1))),
      Selected=c(rep(FALSE,length(phylo_ml_edge_start_values)), identical(ml_selected_name,"topology_search"),
                 seq_along(phylo_ml_edge_start_values)==ml_final_best_i & identical(ml_selected_name,"edge_refined")),
      stringsAsFactors=FALSE
    )
    utils::write.csv(ml_optimisation_diagnostics, file.path(phylo_dir, "ML_staged_optimisation_diagnostics.csv"), row.names=FALSE)
    
    ml_tree <- ape::unroot(ml_fit$tree); ml_tree$node.label <- NULL
    ml_checkpoint_file <- file.path(phylo_dir, paste0("ML_bootstrap_", phylo_ml_model, "_ASC_checkpoint.rds"))
    if (isTRUE(phylo_reset_ml_checkpoint) && file.exists(ml_checkpoint_file)) { message("[ML bootstrap] Removing the existing checkpoint."); file.remove(ml_checkpoint_file) }
    
    ml_bootstrap_trees <- run_ml_bootstrap_progress(ml_fit, phylo_ml_bootstrap_reps, phylo_bootstrap_seed+1L, ml_checkpoint_file,
                                                    phylo_progress_every, phylo_checkpoint_every, phylo_resume_bootstraps, phylo_progress_recent_n)
    message("[ML support] Calculating locus-bootstrap support values."); support_start <- Sys.time()
    ml_supported_tree <- attach_phylo_bootstrap_support(ml_tree, ml_bootstrap_trees)
    message(sprintf("[ML support] Support calculation completed in %s.", format_phylo_duration(difftime(Sys.time(), support_start, units="secs"))))
    
    ape::write.tree(ml_supported_tree, file=file.path(phylo_dir, paste0("Bootstrap_ML_", phylo_ml_model, "_ASC_tree_supported.newick")))
    saveRDS(ml_supported_tree, file=file.path(phylo_dir, paste0("Bootstrap_ML_", phylo_ml_model, "_ASC_tree_supported.rds")))
    saveRDS(ml_fit, file=file.path(phylo_dir, paste0("Bootstrap_ML_", phylo_ml_model, "_ASC_model_fit.rds")))
    if (phylo_save_bootstrap_trees) saveRDS(ml_bootstrap_trees, file=file.path(phylo_dir, "Bootstrap_ML_replicate_trees.rds"))
    message(sprintf("[ML complete] ML analysis completed in %s.", format_phylo_duration(difftime(Sys.time(), ml_stage_start, units="secs"))))
  }
  
  
  #### 5. RECTANGULAR NJ AND ML FIGURES ##################################################
  
  if (!is.null(nj_supported_tree) && !is.null(nj_bootstrap_trees)) {
    message("[NJ outputs] Drawing site- and species-annotated rectangular trees.")
    
    save_supported_phylo_tree(
      nj_supported_tree, nj_bootstrap_trees, site_by_sample, site_phylo_colours, site_phylo_shapes, site_col_name, "species_site",
      "Bootstrap neighbour-joining tree — tips coloured by site", "Bootstrap_NJ_tree_site"
    )
    
    save_supported_phylo_tree(
      nj_supported_tree, nj_bootstrap_trees, species_by_sample, species_phylo_colours, species_phylo_shapes, species_col_name, "species",
      "Bootstrap neighbour-joining tree — tips coloured by species", "Bootstrap_NJ_tree_species"
    )
    
    if (isTRUE(phylo_save_collapsed_species_trees)) {
      save_collapsed_species_phylo_tree(
        nj_supported_tree, nj_bootstrap_trees,
        "Bootstrap neighbour-joining tree — monophyletic species collapsed",
        "Bootstrap_NJ_tree_species_collapsed"
      )
    }
  }
  
  if (!is.null(ml_supported_tree) && !is.null(ml_bootstrap_trees)) {
    message("[ML outputs] Drawing site- and species-annotated rectangular trees.")
    
    save_supported_phylo_tree(
      ml_supported_tree, ml_bootstrap_trees, site_by_sample, site_phylo_colours, site_phylo_shapes, site_col_name, "species_site",
      paste0("Bootstrap ML genotype-state tree (", phylo_ml_model, "+ASC) — tips coloured by site"), "Bootstrap_ML_tree_site"
    )
    
    save_supported_phylo_tree(
      ml_supported_tree, ml_bootstrap_trees, species_by_sample, species_phylo_colours, species_phylo_shapes, species_col_name, "species",
      paste0("Bootstrap ML genotype-state tree (", phylo_ml_model, "+ASC) — tips coloured by species"), "Bootstrap_ML_tree_species"
    )
    
    if (isTRUE(phylo_save_collapsed_species_trees)) {
      save_collapsed_species_phylo_tree(
        ml_supported_tree, ml_bootstrap_trees,
        paste0("Bootstrap ML genotype-state tree (", phylo_ml_model, "+ASC) — monophyletic species collapsed"),
        "Bootstrap_ML_tree_species_collapsed"
      )
    }
  }
  #### 6. TANGLEGRAM FUNCTIONS AND OUTPUTS ################################################
  
  randomly_rotate_tangle_tree <- function(tree, fraction=0.40) {
    internal_nodes <- ape::Ntip(tree)+seq_len(tree$Nnode)
    if (length(internal_nodes)==0L) return(tree)
    selected_nodes <- sample(internal_nodes, min(max(1L,round(length(internal_nodes)*fraction)),length(internal_nodes)), replace=FALSE)
    for (node in selected_nodes) { rotated <- tryCatch(phytools::rotateNodes(tree, nodes=node), error=function(e) NULL); if (!is.null(rotated)) tree <- rotated }
    tree
  }
  
  cophylo_rotation_score <- function(object) {
    if (is.null(object) || !inherits(object,"cophylo")) return(Inf)
    stored <- vapply(object$trees, function(tree) { value <- attr(tree,"minRotate"); if (is.null(value) || length(value)!=1L || !is.finite(value[1])) NA_real_ else as.numeric(value[1]) }, numeric(1))
    if (all(is.finite(stored))) return(sum(stored))
    left_pos <- match(object$assoc[,1], object$trees[[1]]$tip.label); right_pos <- match(object$assoc[,2], object$trees[[2]]$tip.label)
    sum((left_pos-right_pos)^2, na.rm=TRUE)
  }
  
  run_safe_cophylo_rotation <- function(left_tree, right_tree, associations) {
    run_attempt <- function(methods_to_use) {
      tryCatch({
        invisible(utils::capture.output(candidate <- phytools::cophylo(left_tree, right_tree, assoc=associations, rotate=TRUE, methods=methods_to_use,
                                                                       rotate.multi=phylo_tanglegram_rotate_multifurcations, print=FALSE)))
        if (inherits(candidate,"multiCophylo")) {
          candidate_scores <- vapply(candidate, cophylo_rotation_score, numeric(1)); finite_candidates <- which(is.finite(candidate_scores))
          if (length(finite_candidates)==0L) return(list(object=NULL,error="No finite multiCophylo score was produced."))
          candidate <- candidate[[finite_candidates[which.min(candidate_scores[finite_candidates])]]]
        }
        list(object=candidate,error=NA_character_)
      }, error=function(e) list(object=NULL,error=conditionMessage(e)))
    }
    primary <- run_attempt(phylo_tanglegram_methods); if (!is.null(primary$object)) return(primary)
    fallback_method <- as.character(phylo_tanglegram_methods)[1]; fallback <- run_attempt(fallback_method)
    if (!is.null(fallback$object)) { fallback$error <- paste0("Configured methods failed: ",primary$error,". Fallback method '",fallback_method,"' succeeded."); return(fallback) }
    list(object=NULL,error=paste0("Configured rotation failed: ",primary$error," | Fallback rotation failed: ",fallback$error))
  }
  
  optimise_tanglegram_order <- function(left_tree, right_tree, associations, n_starts, random_fraction, seed) {
    set.seed(seed)
    baseline_object <- tryCatch(phytools::cophylo(left_tree, right_tree, assoc=associations, rotate=FALSE), error=function(e) stop("The baseline tanglegram could not be constructed: ",conditionMessage(e),call.=FALSE))
    baseline_score <- cophylo_rotation_score(baseline_object); best_object <- baseline_object; best_score <- baseline_score; selected_start <- 0L
    all_scores <- rep(NA_real_,n_starts); all_errors <- rep(NA_character_,n_starts); all_methods <- rep(NA_character_,n_starts)
    
    for (i in seq_len(n_starts)) {
      message("[Tanglegram rotation] Starting arrangement ",i,"/",n_starts,".")
      left_start <- if (i==1L) left_tree else randomly_rotate_tangle_tree(left_tree,random_fraction)
      right_start <- if (i==1L) right_tree else randomly_rotate_tangle_tree(right_tree,random_fraction)
      attempt <- run_safe_cophylo_rotation(left_start,right_start,associations)
      if (is.null(attempt$object)) { all_errors[i] <- attempt$error; message("[Tanglegram rotation] Arrangement ",i," failed and was skipped: ",attempt$error); next }
      candidate_score <- cophylo_rotation_score(attempt$object); all_scores[i] <- candidate_score
      all_methods[i] <- if (is.na(attempt$error)) paste(phylo_tanglegram_methods,collapse="+") else paste0(as.character(phylo_tanglegram_methods)[1],"_fallback")
      if (!is.na(attempt$error)) all_errors[i] <- attempt$error
      message("[Tanglegram rotation] Arrangement ",i," completed; objective = ",signif(candidate_score,7),".")
      if (is.finite(candidate_score) && candidate_score<best_score) { best_object <- attempt$object; best_score <- candidate_score; selected_start <- i }
    }
    
    if (selected_start==0L) {
      warning("Automatic tanglegram rotation failed or did not improve the baseline. The unrotated arrangement will be used.",call.=FALSE)
      message("[Tanglegram rotation] Using baseline arrangement; objective = ",signif(baseline_score,7),".")
    } else message("[Tanglegram rotation] Selected arrangement ",selected_start," of ",n_starts,"; objective = ",signif(best_score,7),".")
    
    attr(best_object,"rotation_score") <- best_score; attr(best_object,"rotation_scores") <- all_scores; attr(best_object,"rotation_errors") <- all_errors
    attr(best_object,"rotation_methods") <- all_methods; attr(best_object,"selected_start") <- selected_start; attr(best_object,"baseline_score") <- baseline_score
    best_object
  }
  
  prepare_tangle_tree <- function(tree) {
    tree$node.label <- NULL; rooted <- root_phylo_for_display(tree); tree <- ape::collapse.singles(rooted$tree)
    if (phylo_tanglegram_use_branch_lengths) {
      if (is.null(tree$edge.length)) stop("Branch lengths are missing from a tanglegram tree.",call.=FALSE)
      if (any(!is.finite(tree$edge.length))) stop("A tanglegram tree contains non-finite branch lengths.",call.=FALSE)
    } else tree$edge.length <- rep(1,nrow(tree$edge))
    tree <- ape::reorder.phylo(ape::ladderize(tree,right=FALSE),order="cladewise")
    list(tree=tree,method=rooted$method,outgroup_tips=rooted$outgroup_tips)
  }
  
  build_phylo_tanglegram <- function(nj_tree, ml_tree, nj_bootstrap_trees, ml_bootstrap_trees) {
    common_tips <- sort(intersect(nj_tree$tip.label,ml_tree$tip.label)); if (length(common_tips)<4L) stop("Fewer than four shared tips are available.",call.=FALSE)
    nj_rooting <- prepare_tangle_tree(ape::keep.tip(nj_tree,common_tips)); ml_rooting <- prepare_tangle_tree(ape::keep.tip(ml_tree,common_tips))
    if (length(nj_rooting$outgroup_tips)>0L && length(ml_rooting$outgroup_tips)>0L) {
      if (!setequal(nj_rooting$outgroup_tips,ml_rooting$outgroup_tips)) stop("The NJ and ML trees identified different outgroup samples.",call.=FALSE)
      root_method <- nj_rooting$method; message("[Tanglegram rooting] ",root_method," | outgroup samples: ",paste(nj_rooting$outgroup_tips,collapse=", "))
    } else { root_method <- "NJ and ML trees independently midpoint rooted for visualisation only"; message("[Tanglegram rooting] ",root_method) }
    
    associations <- cbind(common_tips,common_tips)
    object <- optimise_tanglegram_order(nj_rooting$tree,ml_rooting$tree,associations,phylo_tanglegram_optimisation_starts,phylo_tanglegram_random_fraction,phylo_tanglegram_rotation_seed)
    nj_bootstrap_common <- prune_bootstrap_trees(nj_bootstrap_trees,common_tips); ml_bootstrap_common <- prune_bootstrap_trees(ml_bootstrap_trees,common_tips)
    
    list(object=object,common_tips=common_tips,root_method=root_method,outgroup_tips=nj_rooting$outgroup_tips,
         left_support=calculate_bootstrap_support(object$trees[[1]],nj_bootstrap_common),right_support=calculate_bootstrap_support(object$trees[[2]],ml_bootstrap_common),
         rotation_score=attr(object,"rotation_score"),rotation_scores=attr(object,"rotation_scores"),rotation_errors=attr(object,"rotation_errors"),
         rotation_methods=attr(object,"rotation_methods"),selected_start=attr(object,"selected_start"),baseline_score=attr(object,"baseline_score"))
  }
  
  draw_phylo_tanglegram <- function(prepared, group_by_sample, colour_palette, shape_palette, group_title, label_mode) {
    label_mode <- match.arg(label_mode,c("species","site","species_site")); object <- prepared$object
    left_samples <- object$trees[[1]]$tip.label; right_samples <- object$trees[[2]]$tip.label
    left_labels <- make_phylo_tip_labels(left_samples,label_mode); right_labels <- make_phylo_tip_labels(right_samples,label_mode)
    link_groups <- clean_phylo_group(unname(group_by_sample[object$assoc[,1]])); link_colours <- unname(colour_palette[link_groups])
    link_colours[is.na(link_colours)] <- "grey60"; link_colours <- grDevices::adjustcolor(link_colours,alpha.f=phylo_tanglegram_link_alpha)
    left_groups <- clean_phylo_group(unname(group_by_sample[left_samples])); right_groups <- clean_phylo_group(unname(group_by_sample[right_samples]))
    left_colours <- unname(colour_palette[left_groups]); right_colours <- unname(colour_palette[right_groups])
    left_shapes <- unname(shape_palette[left_groups]); right_shapes <- unname(shape_palette[right_groups])
    left_colours[is.na(left_colours)] <- "grey50"; right_colours[is.na(right_colours)] <- "grey50"; left_shapes[is.na(left_shapes)] <- 4; right_shapes[is.na(right_shapes)] <- 4
    left_outlines <- ifelse(left_shapes %in% 21:25,"black",left_colours); right_outlines <- ifelse(right_shapes %in% 21:25,"black",right_colours)
    sizes <- get_phylo_sizes(length(prepared$common_tips),tanglegram=TRUE)
    present_groups <- unique(c(left_groups,right_groups)); legend_colours <- unname(colour_palette[present_groups]); legend_shapes <- unname(shape_palette[present_groups])
    legend_colours[is.na(legend_colours)] <- "grey50"; legend_shapes[is.na(legend_shapes)] <- 4; legend_outlines <- ifelse(legend_shapes %in% 21:25,"black",legend_colours)
    legend_layout <- get_phylo_legend_layout(present_groups)
    graphics::layout(matrix(c(1,2),ncol=1),heights=c(1-legend_layout$fraction,legend_layout$fraction))
    graphics::par(mar=c(1.2,1,5.8,1),xpd=NA)
    plot(object,type=c("phylogram","phylogram"),ftype="off",lwd=phylo_tanglegram_tree_width,pts=FALSE,part=0.40,tip.len=0.06,
         tip.lwd=phylo_tanglegram_tip_line_width,tip.lty=phylo_tanglegram_tip_line_type,link.type="straight",link.col=link_colours,
         link.lwd=phylo_tanglegram_link_width,link.lty="solid")
    
    phytools::tiplabels.cophylo(text=left_labels,tip=seq_len(ape::Ntip(object$trees[[1]])),which="left",frame="none",col=left_colours,cex=sizes$label,adj=c(-0.08,0.5))
    phytools::tiplabels.cophylo(text=right_labels,tip=seq_len(ape::Ntip(object$trees[[2]])),which="right",frame="none",col=right_colours,cex=sizes$label,adj=c(1.08,0.5))
    phytools::tiplabels.cophylo(which="left",pch=left_shapes,col=left_outlines,bg=left_colours,cex=sizes$point)
    phytools::tiplabels.cophylo(which="right",pch=right_shapes,col=right_outlines,bg=right_colours,cex=sizes$point)
    
    left_nodes <- which(!is.na(prepared$left_support) & prepared$left_support>=phylo_support_to_plot)
    right_nodes <- which(!is.na(prepared$right_support) & prepared$right_support>=phylo_support_to_plot)
    if (length(left_nodes)>0L) phytools::nodelabels.cophylo(text=format(round(prepared$left_support[left_nodes]),trim=TRUE),node=ape::Ntip(object$trees[[1]])+left_nodes,which="left",frame="none",cex=sizes$support,font=2,col="black")
    if (length(right_nodes)>0L) phytools::nodelabels.cophylo(text=format(round(prepared$right_support[right_nodes]),trim=TRUE),node=ape::Ntip(object$trees[[2]])+right_nodes,which="right",frame="none",cex=sizes$support,font=2,col="black")
    
    rotation_text <- if (is.null(prepared$selected_start) || prepared$selected_start==0L) "Unrotated fallback arrangement" else paste0("Optimised arrangement ",prepared$selected_start)
    graphics::title(main=paste0(species,": NJ versus ML tanglegram — tips coloured by ",group_title,"\nTip labels show ",phylo_tip_label_description(label_mode),
                                "\nLeft: bootstrap NJ | Right: bootstrap ML ",phylo_ml_model,"+ASC\nNode labels are locus-bootstrap support (%) | ",prepared$root_method," | ",rotation_text),cex.main=1.25)
    graphics::par(mar=c(0.2,2,0.2,2),xpd=NA); graphics::plot.new()
    graphics::legend("top",inset=c(0,0.01),legend=present_groups,title=group_title,pch=legend_shapes,col=legend_outlines,pt.bg=legend_colours,
                     pt.cex=1.25,cex=legend_layout$cex,bty="n",ncol=legend_layout$columns,x.intersp=0.55,y.intersp=1.05)
    invisible(object)
  }
  
  save_prepared_tanglegram <- function(prepared, group_by_sample, colour_palette, shape_palette, group_title, label_mode, file_stem) {
    label_mode <- match.arg(label_mode,c("species","site","species_site")); n_tips <- length(prepared$common_tips)
    tip_labels <- make_phylo_tip_labels(prepared$common_tips,label_mode); max_tip_chars <- max(nchar(tip_labels),na.rm=TRUE)
    width_cm <- min(phylo_tanglegram_max_width_cm,max(phylo_tanglegram_min_width_cm,75+max_tip_chars*1.25))
    height_cm <- min(phylo_tanglegram_max_height_cm,max(phylo_tanglegram_min_height_cm,20+n_tips*phylo_tanglegram_height_per_tip_cm))
    message("[Tanglegram output] Saving ",file_stem," at ",round(width_cm,1)," × ",round(height_cm,1)," cm.")
    save_pdf_png(file_stem,width_cm,height_cm,function() draw_phylo_tanglegram(prepared,group_by_sample,colour_palette,shape_palette,group_title,label_mode))
  }
  
  if (phylo_run_tanglegrams) {
    if (is.null(nj_supported_tree) || is.null(ml_supported_tree) || is.null(nj_bootstrap_trees) || is.null(ml_bootstrap_trees)) {
      message("[Tanglegram] Not generated because completed NJ and ML analyses and bootstrap trees are required.")
    } else {
      prepared_tanglegram <- build_phylo_tanglegram(nj_supported_tree,ml_supported_tree,nj_bootstrap_trees,ml_bootstrap_trees)
      saveRDS(prepared_tanglegram,file=file.path(phylo_dir,"Tanglegram_NJ_vs_ML_optimised_object.rds"))
      
      rotation_table <- data.frame(Starting_arrangement=seq_along(prepared_tanglegram$rotation_scores),Rotation_objective=prepared_tanglegram$rotation_scores,
                                   Rotation_method=prepared_tanglegram$rotation_methods,Rotation_error=prepared_tanglegram$rotation_errors,Selected=FALSE,stringsAsFactors=FALSE)
      if (!is.null(prepared_tanglegram$selected_start) && prepared_tanglegram$selected_start>0L) rotation_table$Selected[prepared_tanglegram$selected_start] <- TRUE
      utils::write.csv(rotation_table,file.path(phylo_dir,"Tanglegram_NJ_vs_ML_rotation_scores.csv"),row.names=FALSE)
      utils::write.csv(data.frame(Baseline_rotation_objective=prepared_tanglegram$baseline_score,Selected_rotation_objective=prepared_tanglegram$rotation_score,
                                  Selected_start=prepared_tanglegram$selected_start,Used_unrotated_fallback=prepared_tanglegram$selected_start==0L),
                       file.path(phylo_dir,"Tanglegram_NJ_vs_ML_rotation_summary.csv"),row.names=FALSE)
      
      save_prepared_tanglegram(prepared_tanglegram,site_by_sample,site_phylo_colours,site_phylo_shapes,site_col_name,"species_site","Tanglegram_NJ_vs_ML_site")
      save_prepared_tanglegram(prepared_tanglegram,species_by_sample,species_phylo_colours,species_phylo_shapes,species_col_name,"species","Tanglegram_NJ_vs_ML_species")
    }
  }
  
  
  #### 7. MISSING-DATA ANNOTATED TREES AND TANGLEGRAM ############################
  
  if (phylo_generate_missing_data_plots) {
    phylo_missing_count <- setNames(rowSums(is.na(phylo_gt)),rownames(phylo_gt))
    phylo_missing_percent <- setNames(100*rowMeans(is.na(phylo_gt)),rownames(phylo_gt))
    phylo_missing_actual_range <- range(phylo_missing_percent,finite=TRUE); phylo_missing_scale_range <- phylo_missing_actual_range
    if (!is.finite(diff(phylo_missing_scale_range)) || diff(phylo_missing_scale_range)==0) phylo_missing_scale_range <- phylo_missing_scale_range+c(-0.5,0.5)
    phylo_missing_palette <- grDevices::hcl.colors(101,palette=phylo_missing_palette_name,rev=isTRUE(phylo_missing_high_values_dark))
    
    phylo_missing_colours <- function(values) {
      scaled <- (values-phylo_missing_scale_range[1])/diff(phylo_missing_scale_range); scaled <- pmin(1,pmax(0,scaled))
      phylo_missing_palette[1L+round(scaled*(length(phylo_missing_palette)-1L))]
    }
    
    phylo_missing_labels <- function(samples) {
      base_labels <- make_phylo_tip_labels(samples,"species_site")
      paste0(base_labels," [",formatC(phylo_missing_percent[samples],format="f",digits=phylo_missing_label_digits),"% missing]")
    }
    
    draw_phylo_missing_legend <- function() {
      graphics::par(mar=c(0.4,2,0.1,2),xpd=NA); graphics::plot.new(); graphics::plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
      x_positions <- seq(0.15,0.85,length.out=length(phylo_missing_palette)+1L)
      graphics::rect(x_positions[-length(x_positions)],0.42,x_positions[-1L],0.62,col=phylo_missing_palette,border=NA); graphics::rect(0.15,0.42,0.85,0.62,border="black")
      legend_ticks <- if (diff(phylo_missing_actual_range)==0) phylo_missing_actual_range[1] else pretty(phylo_missing_actual_range,n=5)
      legend_ticks <- legend_ticks[legend_ticks>=phylo_missing_actual_range[1] & legend_ticks<=phylo_missing_actual_range[2]]
      tick_positions <- 0.15+(legend_ticks-phylo_missing_scale_range[1])/diff(phylo_missing_scale_range)*0.70
      graphics::segments(tick_positions,0.37,tick_positions,0.42); graphics::text(tick_positions,0.29,labels=formatC(legend_ticks,format="f",digits=phylo_missing_label_digits),cex=0.78)
      graphics::text(0.50,0.81,"Missing genotype data across analysed loci (%)",font=2,cex=0.90)
    }
    
    draw_supported_phylo_missing <- function(tree,bootstrap_trees,method_title) {
      display_rooting <- prepare_supported_display_tree(tree); display_tree <- display_rooting$tree
      support <- calculate_bootstrap_support(display_tree,bootstrap_trees); sample_names <- display_tree$tip.label
      tip_colours <- phylo_missing_colours(unname(phylo_missing_percent[sample_names])); label_colours <- if (phylo_tree_colour_sample_labels) tip_colours else rep("black",length(tip_colours))
      plot_tree <- display_tree; plot_tree$tip.label <- phylo_missing_labels(sample_names); sizes <- get_phylo_sizes(ape::Ntip(plot_tree))
      tree_span <- if (!phylo_tree_show_branch_lengths || is.null(plot_tree$edge.length)) 1 else max(ape::node.depth.edgelength(plot_tree),na.rm=TRUE)
      if (!is.finite(tree_span) || tree_span<=0) tree_span <- 1
      graphics::layout(matrix(c(1,2),ncol=1),heights=c(0.94,0.06)); graphics::par(mar=c(2.2,2,5.0,3),xpd=NA)
      ape::plot.phylo(plot_tree,type="phylogram",direction="rightwards",use.edge.length=phylo_tree_show_branch_lengths,show.tip.label=TRUE,tip.color=label_colours,
                      cex=sizes$label,edge.width=0.70,label.offset=tree_span*0.008,align.tip.label=FALSE,no.margin=FALSE)
      ape::tiplabels(pch=21,col="black",bg=tip_colours,cex=sizes$point,lwd=0.7)
      supported_nodes <- which(!is.na(support) & support>=phylo_support_to_plot)
      if (length(supported_nodes)>0L) ape::nodelabels(text=format(round(support[supported_nodes]),trim=TRUE),node=ape::Ntip(plot_tree)+supported_nodes,frame="none",cex=sizes$support,font=2,col="black")
      graphics::title(main=paste0(species,": ",method_title,"\nTip labels show species, site, accession and missing data",
                                  "\nNode labels are locus-bootstrap support (%) | ",display_rooting$method),cex.main=1.25)
      draw_phylo_missing_legend()
    }
    
    save_supported_phylo_missing <- function(tree,bootstrap_trees,method_title,file_stem) {
      labels <- phylo_missing_labels(tree$tip.label); n_tips <- ape::Ntip(tree); max_label_chars <- max(nchar(labels),na.rm=TRUE)
      width_cm <- min(phylo_tree_max_width_cm,max(phylo_tree_min_width_cm,50+max_label_chars*1.00))
      height_cm <- min(phylo_tree_max_height_cm,max(phylo_tree_min_height_cm,20+n_tips*phylo_tree_height_per_tip_cm))
      save_pdf_png(file_stem,width_cm,height_cm,function() draw_supported_phylo_missing(tree,bootstrap_trees,method_title))
    }
    
    draw_phylo_tanglegram_missing <- function(prepared) {
      object <- prepared$object; link_colours <- phylo_missing_colours(phylo_missing_percent[object$assoc[,1]])
      link_colours <- grDevices::adjustcolor(link_colours,alpha.f=phylo_tanglegram_link_alpha)
      left_samples <- object$trees[[1]]$tip.label; right_samples <- object$trees[[2]]$tip.label
      left_colours <- phylo_missing_colours(phylo_missing_percent[left_samples]); right_colours <- phylo_missing_colours(phylo_missing_percent[right_samples])
      sizes <- get_phylo_sizes(length(prepared$common_tips),tanglegram=TRUE)
      graphics::layout(matrix(c(1,2),ncol=1),heights=c(0.94,0.06)); graphics::par(mar=c(1.2,1,5.5,1),xpd=NA)
      plot(object,type=c("phylogram","phylogram"),ftype="off",lwd=phylo_tanglegram_tree_width,pts=FALSE,part=0.40,tip.len=0.06,
           tip.lwd=phylo_tanglegram_tip_line_width,tip.lty=phylo_tanglegram_tip_line_type,link.type="straight",link.col=link_colours,
           link.lwd=phylo_tanglegram_link_width,link.lty="solid")
      phytools::tiplabels.cophylo(text=phylo_missing_labels(left_samples),tip=seq_len(ape::Ntip(object$trees[[1]])),which="left",frame="none",col=left_colours,cex=sizes$label,adj=c(-0.08,0.5))
      phytools::tiplabels.cophylo(text=phylo_missing_labels(right_samples),tip=seq_len(ape::Ntip(object$trees[[2]])),which="right",frame="none",col=right_colours,cex=sizes$label,adj=c(1.08,0.5))
      phytools::tiplabels.cophylo(which="left",pch=21,col="black",bg=left_colours,cex=sizes$point); phytools::tiplabels.cophylo(which="right",pch=21,col="black",bg=right_colours,cex=sizes$point)
      left_nodes <- which(!is.na(prepared$left_support) & prepared$left_support>=phylo_support_to_plot); right_nodes <- which(!is.na(prepared$right_support) & prepared$right_support>=phylo_support_to_plot)
      if (length(left_nodes)>0L) phytools::nodelabels.cophylo(text=format(round(prepared$left_support[left_nodes]),trim=TRUE),node=ape::Ntip(object$trees[[1]])+left_nodes,which="left",frame="none",cex=sizes$support,font=2,col="black")
      if (length(right_nodes)>0L) phytools::nodelabels.cophylo(text=format(round(prepared$right_support[right_nodes]),trim=TRUE),node=ape::Ntip(object$trees[[2]])+right_nodes,which="right",frame="none",cex=sizes$support,font=2,col="black")
      graphics::title(main=paste0(species,": NJ versus ML tanglegram — missing genotype data\nLeft: bootstrap NJ | Right: bootstrap ML ",phylo_ml_model,
                                  "+ASC\nTip labels show species, site, accession and missing data | ",prepared$root_method),cex.main=1.25)
      draw_phylo_missing_legend(); invisible(object)
    }
    
    save_phylo_tanglegram_missing <- function(prepared,file_stem) {
      labels <- phylo_missing_labels(prepared$common_tips); n_tips <- length(prepared$common_tips); max_label_chars <- max(nchar(labels),na.rm=TRUE)
      width_cm <- min(phylo_tanglegram_max_width_cm,max(phylo_tanglegram_min_width_cm,80+max_label_chars*1.25))
      height_cm <- min(phylo_tanglegram_max_height_cm,max(phylo_tanglegram_min_height_cm,20+n_tips*phylo_tanglegram_height_per_tip_cm))
      save_pdf_png(file_stem,width_cm,height_cm,function() draw_phylo_tanglegram_missing(prepared))
    }
    
    phylo_missing_table <- data.frame(sample=rownames(phylo_gt),missing_genotypes=unname(phylo_missing_count[rownames(phylo_gt)]),analysed_loci=ncol(phylo_gt),
                                      missing_percent=unname(phylo_missing_percent[rownames(phylo_gt)]),site=unname(site_by_sample[rownames(phylo_gt)]),
                                      species=unname(species_by_sample[rownames(phylo_gt)]),stringsAsFactors=FALSE)
    phylo_missing_table <- phylo_missing_table[order(phylo_missing_table$missing_percent,decreasing=TRUE),]
    utils::write.csv(phylo_missing_table,file.path(phylo_dir,"Phylogenetic_tip_missing_data.csv"),row.names=FALSE)
    
    if (!is.null(nj_supported_tree) && !is.null(nj_bootstrap_trees)) {
      message("[Missing-data outputs] Drawing missing-data annotated NJ tree.")
      save_supported_phylo_missing(nj_supported_tree,nj_bootstrap_trees,"Bootstrap neighbour-joining tree — tips by missing data","Bootstrap_NJ_tree_missing_data")
    }
    if (!is.null(ml_supported_tree) && !is.null(ml_bootstrap_trees)) {
      message("[Missing-data outputs] Drawing missing-data annotated ML tree.")
      save_supported_phylo_missing(ml_supported_tree,ml_bootstrap_trees,paste0("Bootstrap ML genotype-state tree (",phylo_ml_model,"+ASC) — tips by missing data"),"Bootstrap_ML_tree_missing_data")
    }
    if (!is.null(prepared_tanglegram)) {
      message("[Missing-data outputs] Drawing missing-data annotated NJ-versus-ML tanglegram.")
      save_phylo_tanglegram_missing(prepared_tanglegram,"Tanglegram_NJ_vs_ML_missing_data")
    }
    message("[Missing-data outputs] Missing-data phylogenies and sample table written to: ",phylo_dir)
  }
  
  
  #### 8. SUMMARY TABLES AND SESSION INFORMATION #########################################
  
  support_tables <- list()
  if (!is.null(nj_supported_tree)) support_tables[["NJ"]] <- extract_support_table(nj_supported_tree,"Neighbour joining")
  if (!is.null(ml_supported_tree)) support_tables[["ML"]] <- extract_support_table(ml_supported_tree,paste0("ML ",phylo_ml_model,"+ASC"))
  if (length(support_tables)>0L) utils::write.csv(do.call(rbind,support_tables),file.path(phylo_dir,"Phylogenetic_node_support_values.csv"),row.names=FALSE)
  
  comparison_summary <- data.frame(
    Samples=nrow(phylo_gt),Polymorphic_loci=ncol(phylo_gt),
    NJ_bootstrap_replicates=if (!is.null(nj_bootstrap_trees)) length(nj_bootstrap_trees) else NA_integer_,
    ML_bootstrap_replicates=if (!is.null(ml_bootstrap_trees)) length(ml_bootstrap_trees) else NA_integer_,
    ML_model=if (!is.null(ml_supported_tree)) paste0(phylo_ml_model,"+ASC") else NA_character_,
    ML_rearrangement=if (!is.null(ml_supported_tree)) phylo_ml_rearrangement else NA_character_,
    ML_log_likelihood=if (!is.null(ml_fit)) as.numeric(ml_fit$logLik) else NA_real_,
    ML_maximum_edge=if (!is.null(ml_fit)) max(ml_fit$tree$edge.length) else NA_real_,
    ML_total_tree_length=if (!is.null(ml_fit)) sum(ml_fit$tree$edge.length) else NA_real_,
    Rooting=if (length(configured_outgroup_taxa)>0L) paste0("Outgroup: ",paste(configured_outgroup_taxa,collapse=", ")) else "Midpoint display rooting",
    Tree_branch_display=if (phylo_tree_show_branch_lengths) "Phylogram: analytical branch lengths shown" else "Cladogram: equal-depth topology display",
    Site_figure_tip_labels="species | site | sample accession",Species_figure_tip_labels="species | sample accession",Replot_only=phylo_replot_only,
    stringsAsFactors=FALSE
  )
  
  if (!is.null(nj_supported_tree) && !is.null(ml_supported_tree)) comparison_summary$Normalised_Robinson_Foulds_distance <- tryCatch(
    as.numeric(phangorn::RF.dist(nj_supported_tree,ml_supported_tree,normalize=TRUE)),error=function(e) NA_real_)
  
  utils::write.csv(comparison_summary,file.path(phylo_dir,"Phylogenetic_analysis_summary.csv"),row.names=FALSE)
  utils::write.csv(data.frame(sample=phylo_meta$.phylo_sample,site=unname(site_by_sample),species=unname(species_by_sample),stringsAsFactors=FALSE),
                   file.path(phylo_dir,"Phylogenetic_tip_metadata.csv"),row.names=FALSE)
  capture.output(utils::sessionInfo(),file=file.path(phylo_dir,"Phylogenetic_analysis_sessionInfo.txt"))
  message("Bootstrap-supported phylogenetic outputs written to: ",phylo_dir)
}
