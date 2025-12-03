report.dart.qc.stats.RD <- function (dart_data, basedir, species, dataset, threshold_missing_loci = 0.5) 
{
  name <- dataset
  genotypes <- dart_data$gt
  treatment <- dart_data$treatment
  locus_repro <- dart_data$locus_repro
  num_samples <- nrow(genotypes)
  num_loci <- ncol(genotypes)
  num_clones <- length(unique(dart_data$locus_names))
  qualstat_dir <- paste(basedir, "/qual_stat/", 
                        sep = "")
  
  if (!dir.exists(qualstat_dir)) {
    cat("  QC directory: ", qualstat_dir, " does not exist and is being created. \n")
    dir.create(qualstat_dir)
  }
  
  qc_dir <- paste(basedir, "/qual_stat/", treatment, 
                  sep = "")
  
  if (!dir.exists(qc_dir)) {
    cat("  QC directory: ", qc_dir, " does not exist and is being created. \n")
    dir.create(qc_dir)
  }
  
  else {
    cat("  QC directory: ", qc_dir, " already exists, content will be overwritten. \n")
  }
  if (exists("genotypes") & num_samples >= 1 & num_loci >= 
      1) {
    cat("  QC calculations starting \n")
  }
  else {
    cat(" Fatal Error: could not identify a genotype table...\n")
    stop()
  }
  missing <- is.na(dart_data$gt)
  count_of_missing_by_locus <- rowSums(missing)
  count_of_missing_by_sample <- colSums(missing)
  missing_lines <- c(paste("Samples missing data for ", threshold_missing_loci, 
                           " of loci"))
  ind_samples_high_missing <- which((count_of_missing_by_locus/num_loci) > 
                                      threshold_missing_loci)
  if (length(ind_samples_high_missing) >= 1) {
    for (i in 1:length(ind_samples_high_missing)) {
      sample_name_i <- rownames(genotypes)[ind_samples_high_missing[i]]
      prop_missing_i <- (count_of_missing_by_locus/num_loci)[ind_samples_high_missing[i]]
      missing_lines <- c(missing_lines, paste(sample_name_i, 
                                              prop_missing_i, sep = ": "))
    }
    
samples_missing_file <- paste(basedir, "/qual_stat/", treatment, "/samples_exceeding_missing_threshold.txt", sep = "")

    write.table(rownames(genotypes)[ind_samples_high_missing], 
                file = samples_missing_file, sep = ",", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
    report_missing <- samples_missing_file
  }
  else {
    missing_lines <- c(missing_lines, "None")
    report_missing <- "None"
  }
  AF.summary <- calculate.AF.summary(dart_data)
  pdf(file = paste(basedir, species, "/qual_stat/", treatment, 
                   "/missing_data_by_sample.pdf", sep = ""))
  hist(count_of_missing_by_sample/num_samples, xlab = "proportion of samples that are missing", 
       ylab = "frequency", main = "")
  dev.off()
  pdf(file = paste(basedir, species, "/qual_stat/", treatment, 
                   "/missing_data_by_locus.pdf", sep = ""))
  hist(count_of_missing_by_locus/num_loci, 200, xlab = "proportion of loci that are missing", 
       ylab = "frequency", main = "")
  dev.off()
  pdf(file = paste(basedir, species, "/qual_stat/", treatment, 
                   "/reprodicibility_of_loci.pdf", sep = ""))
  hist(locus_repro, 200, xlab = "locus reproducibility score", 
       ylab = "frequency", main = "")
  dev.off()
  pdf(file = paste(basedir, species, "/qual_stat/", treatment, 
                   "/AF_versus_He.pdf", sep = ""))
  plot(AF.summary$P, AF.summary$H, xlab = "Allele Frequency", 
       ylab = "Heterozygosity", main = "")
  dev.off()
  date_line <- date()
  species_line <- paste("Species: ", species, sep = "")
  name_line <- paste("Name: ", name, sep = "")
  snps_line <- paste("Number of SNPs: ", num_loci, sep = "")
  clones_line <- paste("Number of clones: ", num_clones, sep = "")
  samples_line <- paste("Number of samples: ", num_samples, 
                        sep = "")
  report_file <- file(paste(basedir, species, "/qual_stat/", 
                            treatment, "/quality_report.txt", sep = ""))
  report <- c(date_line, species_line, name_line, snps_line, 
              clones_line, samples_line, "\n", missing_lines, "\n")
  writeLines(report, report_file)
  close(report_file)
  cat(paste("  QC report in ", basedir, species, "/qual_stat/", 
            treatment, "/\n", sep = ""))
  QC_report = list(treatment = treatment, report_missing = report_missing, 
                   threshold_missing_loci = threshold_missing_loci)
  return(QC_report)
}

