#functions for running allele proportion scripts - R. Dimon 2024


dist_row_sums <- function(xy, n){
  subset <- xy
  alldist <- as.matrix(dist(subset))
  while (nrow(subset) > n) {
    cdists = rowSums(alldist)
    closest <- which(cdists == min(cdists))[1]
    subset <- subset[-closest,]
    alldist <- alldist[-closest,-closest]
  }
  return(subset)
}







Common_Allele_Prop_Random_Sites <- function(dms = dms, gt_sw_comp=gt_sw_comp, analysis=analysis, NumSteps = NumSteps, i_sw_common=i_sw_common, samplethreshold = samplethreshold, sitethreshold = sitethreshold){
  
  # dms = data
  # NumSteps = Number of Loops
  # i_sw_common = common allele dataset (minor allele frequencies greater than MAF threshold)
  # samplethreshold = number of individuals to sample per site
  # sitethreshold= max number of sites to obtain sampling combination for (high numnber of sites is redundant as its becomes very close to 100% common alleles)

    sampling<-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("sites", "samples", "total"))
    for (c in 1:length(unique(dms$meta$analyses[,analysis]))) {
      for (ind in 1:samplethreshold) {
        sampling[nrow(sampling) + 1,] = list(c,ind,c*ind)
      }
    }
    sampling <- sampling[which(sampling$samples==samplethreshold & sampling$sites<=sitethreshold | sampling$samples==(samplethreshold-1) & sampling$sites<=sitethreshold | sampling$samples==(samplethreshold-2) & sampling$sites<=sitethreshold),]
    sampling<-sampling[order(sampling$total),]
    dms_meta <- cbind.data.frame(sample=dms$meta$sample_names, site=dms$meta$analyses[,analysis],lat=dms$meta$lat, long=dms$meta$long)
    colnames(dms_meta) <- c("sample","site","lat","long")
    unique_sites <- unique(na.omit(dms_meta$site)) #maximum number of sites
    gvals <- list() #accumulate data to be saved #samples & sites to choose, number of common alleles
    #insert loop here to run through table listing parameters in bottom two lines#
    OGM_DF<-setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("m","t_num_indv", "n_sites_sel", "n_indiv_sel", "rand_sites_sel", "rand_indiv_sel", "jvals", "Aprop"))
    set.seed(12345)
    
    for (n in 1:length(sampling$sites)) {
      n_sites_sel <- sampling$sites[n]
      n_indiv_sel <- sampling$samples[n]
      cat("running...", n,"... ",n_sites_sel," sites...", n_indiv_sel," samples...\n")
      
      t_num_indv <- n_sites_sel*n_indiv_sel
      for (m in 1:NumSteps) {
        rand_sites_sel <- unique_sites[sample(1:length(unique_sites))[1:n_sites_sel]]; #print(rand_sites_sel)
        rand_indiv_sel <- c()
        for (s in 1:n_sites_sel) {
          rand_indiv <- as.character(dms_meta$sample[which(dms_meta$site == rand_sites_sel[s])][sample(1:length(which(dms_meta$site == rand_sites_sel[s])))[1:n_indiv_sel]])
          rand_indiv_sel <- c(rand_indiv_sel,rand_indiv) ; #print(rand_indiv_sel)
        }
        fixed_indi <- which(dms$sample_names %in% rand_indiv_sel)
        ran_vec <- rep(0, nrow(gt_sw_comp))
        ran_vec[fixed_indi] <- 1
        common_alleles  <- common_allele_count(gt_sw_comp, ran_vec)
        jvals <- length(intersect( which(common_alleles[[2]] > 0), i_sw_common))
        Aprop<-jvals/length(i_sw_common)
        gvals[[m]] <- data.frame(m,t_num_indv,n_sites_sel,n_indiv_sel,rand_sites_sel,rand_indiv_sel,jvals, Aprop)
      }
      big_data2 = do.call(rbind, gvals)
      OGM_DF<-rbind(OGM_DF, big_data2)
    }
    clean_data <- OGM_DF %>% group_by(t_num_indv,n_sites_sel,m) %>% summarize(Allele=mean(Aprop),Allelen=mean(jvals), indv_p_site = n_indiv_sel) %>% data.frame
    clean_data$t_num_indv<-as.numeric(clean_data$t_num_indv)
    # clean_data$t_num_indv <- paste("total=",clean_data$t_num_indv)
    # clean_data$n_sites_sel<-paste(clean_data$n_sites_sel,"sites")
    # clean_data$indv_p_site<-paste(clean_data$indv_p_site,"indv")
    
    
    
    return(clean_data)
    
   
}






Common_Allele_Prop_Fixed_Sites <- function(dms = dms, gt_sw_comp=gt_sw_comp, analysis=analysis, i_sw_common=i_sw_common, samplethreshold = samplethreshold, fixedsites = fixedsites){

  # dms = data
  # i_sw_common = common allele dataset (minor allele frequencies greater than MAF threshold)
  # samplethreshold = number of individuals to sample per site
  # fixedsites = unique site names you want to obtain common allele prop
  

    dms_meta <- cbind.data.frame(sample=dms$meta$sample_names, site=dms$meta$analyses[,analysis],lat=dms$meta$lat, long=dms$meta$long)
    colnames(dms_meta) <- c("sample","site","lat","long")
    
    
    
    rand_indiv_sel <- c()
    
    for (s in 1:length(fixedsites)) {
      rand_indiv <- as.character(dms_meta$sample[which(dms_meta$site == fixedsites[s])][sample(1:length(which(dms_meta$site == fixedsites[s])))[1:samplethreshold]])
      rand_indiv_sel <- c(rand_indiv_sel,rand_indiv) ; #print(rand_indiv_sel)
    }
    
    fixed_indi <- which(dms$sample_names %in% rand_indiv_sel)
    ran_vec <- rep(0, nrow(gt_sw_comp))
    ran_vec[fixed_indi] <- 1
    common_alleles  <- common_allele_count(gt_sw_comp, ran_vec)
    jvals <- length(intersect( which(common_alleles[[2]] > 0), i_sw_common))
    Aprop<-jvals/length(i_sw_common)
    
    return(Aprop)

}






Common_Allele_Prop_Median_Fixed_Sites <- function(dms = dms, gt_sw_comp=gt_sw_comp, analysis=analysis, i_sw_common=i_sw_common, samplethreshold = samplethreshold, sitethreshold = sitethreshold){
  
  #identifies sites with the greatest median distance from each other and calculates common alleles
  
  # dms = data
  # i_sw_common = common allele dataset (minor allele frequencies greater than MAF threshold)
  # samplethreshold = number of individuals to sample per site
  # sitethreshold = number of sites to sample combinations up to
  
  require(matrixStats)
  require(dplyr)
  
  dms_meta <- cbind.data.frame(sample=dms$meta$sample_names, site=dms$meta$analyses[,analysis],lat=dms$meta$lat, long=dms$meta$long)
  colnames(dms_meta) <- c("sample","site","lat","long")
  avgsites <- dms_meta %>% group_by(site) %>% summarize(alat=mean(lat),along=mean(long)) %>% data.frame
  
  
  
  sampling<-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("sites", "samples", "total"))
  for (c in 1:length(unique(dms$meta$analyses[,analysis]))) {
    for (ind in 1:samplethreshold) {
      sampling[nrow(sampling) + 1,] = list(c,ind,c*ind)
    }
  }
  sampling <- sampling[which(sampling$samples==samplethreshold & sampling$sites<=sitethreshold | sampling$samples==(samplethreshold-1) & sampling$sites<=sitethreshold),]
  sampling<-sampling[order(sampling$total),]
  
  siteCombos <- sampling$sites
  
  
  
  dist_row_medians <- function(xy, n){
    subset <- xy
    alldist <- as.matrix(dist(subset[-1]))
    while (nrow(subset) > n) {
      cdists = rowMedians(alldist)
      closest <- which(cdists == min(cdists))[1]
      subset <- subset[-closest,]
      alldist <- alldist[-closest,-closest]
    }
    return(subset)
  }
  
  MedFixProp_df <- data.frame(matrix(ncol=6, nrow=length(siteCombos)))
  colnames(MedFixProp_df) <- c("t_num_indv", "n_sites_sel", "Allele", "Allelen", "samples", "SiteNames") 
  
  MedFixProp_df$t_num_indv <- samplethreshold
  MedFixProp_df$n_sites_sel <- siteCombos
  MedFixProp_df$samples <- samplethreshold*siteCombos
  
  
  for (n in siteCombos){
    
  numb <- siteCombos[n]
  mediansites = dist_row_medians(avgsites,numb)
  fixedsites <- mediansites$site
  
  rand_indiv_sel <- c()
  
  for (s in 1:length(fixedsites)) {
    rand_indiv <- as.character(dms_meta$sample[which(dms_meta$site == fixedsites[s])][sample(1:length(which(dms_meta$site == fixedsites[s])))[1:samplethreshold]])
    rand_indiv_sel <- c(rand_indiv_sel,rand_indiv) ; #print(rand_indiv_sel)
  }
  
  fixed_indi <- which(dms$sample_names %in% rand_indiv_sel)
  ran_vec <- rep(0, nrow(gt_sw_comp))
  ran_vec[fixed_indi] <- 1
  common_alleles  <- common_allele_count(gt_sw_comp, ran_vec)
  jvals <- length(intersect( which(common_alleles[[2]] > 0), i_sw_common))
  Aprop<-jvals/length(i_sw_common)
  
  MedFixProp_df$Allelen[n] <- jvals
  MedFixProp_df$Allele[n] <- Aprop
  ufixedsites <- unlist(fixedsites)
  MedFixProp_df$SiteNames[n] <- paste(ufixedsites, collapse=" ")
  
  }
  
  
  return(MedFixProp_df)
  
}




