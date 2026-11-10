report.dart.qc.stats.RD <- function(dart_data, basedir, species, dataset, threshold_missing_loci = 0.5) {
  
  # Quality-control report for a DArT object, written to <basedir>/qual_stat/<treatment>/.
  #
  # Files written:
  #   missing_data_by_sample.pdf  proportion of loci missing in each sample
  #   missing_data_by_locus.pdf   proportion of samples missing at each locus
  #   reprodicibility_of_loci.pdf DArT reproducibility score of each locus
  #   AF_versus_Ho.pdf            allele frequency vs observed heterozygosity, with
  #                               the Hardy-Weinberg expectation 2p(1-p) and the
  #                               maximum possible heterozygosity 2*min(p, 1-p)
  #   quality_report.txt          counts and samples above threshold_missing_loci
  #   samples_exceeding_missing_threshold.txt (only when any sample exceeds it)
  #
  # Rows of dart_data$gt are samples and columns are loci.
  
  name <- dataset
  genotypes <- dart_data$gt
  treatment <- dart_data$treatment
  locus_repro <- dart_data$locus_repro
  
  if (is.null(genotypes) || nrow(genotypes) < 1 || ncol(genotypes) < 1) {
    stop("report.dart.qc.stats.RD: could not identify a genotype table in dart_data$gt.", call. = FALSE)
  }
  
  num_samples <- nrow(genotypes)
  num_loci <- ncol(genotypes)
  num_clone_ids <- length(unique(dart_data$locus_names))   # distinct DArT clone IDs (sequence tags), not clonal individuals
  
  qualstat_dir <- file.path(basedir, "qual_stat")
  qc_dir <- file.path(qualstat_dir, treatment)
  if (!dir.exists(qualstat_dir)) {
    cat("  QC directory: ", qualstat_dir, " does not exist and is being created. \n")
  }
  if (!dir.exists(qc_dir)) {
    cat("  QC directory: ", qc_dir, " does not exist and is being created. \n")
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    cat("  QC directory: ", qc_dir, " already exists, content will be overwritten. \n")
  }
  cat("  QC calculations starting \n")
  
  # Missing data. Rows are samples and columns are loci.
  missing <- is.na(genotypes)
  prop_missing_per_sample <- rowSums(missing) / num_loci
  prop_missing_per_locus <- colSums(missing) / num_samples
  
  # Samples above the missingness threshold.
  missing_lines <- paste("Samples missing data for more than", threshold_missing_loci, "of loci")
  ind_samples_high_missing <- which(prop_missing_per_sample > threshold_missing_loci)
  if (length(ind_samples_high_missing) >= 1) {
    missing_lines <- c(
      missing_lines,
      paste(rownames(genotypes)[ind_samples_high_missing],
            round(prop_missing_per_sample[ind_samples_high_missing], 4), sep = ": ")
    )
    samples_missing_file <- file.path(qc_dir, "samples_exceeding_missing_threshold.txt")
    write.table(rownames(genotypes)[ind_samples_high_missing], file = samples_missing_file, sep = ",",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    report_missing <- samples_missing_file
  } else {
    missing_lines <- c(missing_lines, "None")
    report_missing <- "None"
  }
  
  AF.summary <- calculate.AF.summary(dart_data)
  
  # Proportion of loci missing in each sample (one value per sample).
  pdf(file = file.path(qc_dir, "missing_data_by_sample.pdf"))
  hist(prop_missing_per_sample, breaks = 50,
       xlab = "Proportion of loci missing in each sample", ylab = "Number of samples",
       main = paste0("Missing data by sample (n = ", num_samples, " samples)"))
  abline(v = threshold_missing_loci, lty = 2, col = "red")
  dev.off()
  
  # Proportion of samples missing at each locus (one value per locus).
  pdf(file = file.path(qc_dir, "missing_data_by_locus.pdf"))
  hist(prop_missing_per_locus, breaks = 200,
       xlab = "Proportion of samples missing at each locus", ylab = "Number of loci",
       main = paste0("Missing data by locus (n = ", num_loci, " loci)"))
  dev.off()
  
  pdf(file = file.path(qc_dir, "reprodicibility_of_loci.pdf"))
  hist(locus_repro, 200, xlab = "Locus reproducibility score", ylab = "Number of loci", main = "")
  dev.off()
  
  # Allele frequency vs observed heterozygosity (H = proportion of heterozygous calls).
  pdf(file = file.path(qc_dir, "AF_versus_Ho.pdf"))
  plot(AF.summary$P, AF.summary$H, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Allele frequency", ylab = "Observed heterozygosity (Ho)",
       main = paste0(num_samples, " samples, ", num_loci, " loci"))
  p_grid <- seq(0, 1, length.out = 201)
  lines(p_grid, 2 * p_grid * (1 - p_grid), col = "blue", lwd = 2)
  lines(p_grid, 2 * pmin(p_grid, 1 - p_grid), col = "red", lwd = 1.5, lty = 2)
  legend("topright", bty = "n", cex = 0.8, lwd = c(2, 1.5), lty = c(1, 2), col = c("blue", "red"),
         legend = c("Hardy-Weinberg expectation 2p(1-p)", "Maximum possible Ho 2*min(p, 1-p)"))
  dev.off()
  
  report <- c(
    date(),
    paste0("Species: ", species),
    paste0("Name: ", name),
    paste0("Number of SNPs: ", num_loci),
    paste0("Number of DArT clone IDs (loci): ", num_clone_ids),
    paste0("Number of samples: ", num_samples),
    paste0("Mean proportion of loci missing per sample: ", round(mean(prop_missing_per_sample), 4)),
    paste0("Mean proportion of samples missing per locus: ", round(mean(prop_missing_per_locus), 4)),
    "\n",
    missing_lines,
    "\n"
  )
  report_file <- file(file.path(qc_dir, "quality_report.txt"))
  writeLines(report, report_file)
  close(report_file)
  cat(paste0("  QC report in ", qc_dir, "/\n"))
  
  QC_report <- list(treatment = treatment, report_missing = report_missing, threshold_missing_loci = threshold_missing_loci)
  return(QC_report)
}
