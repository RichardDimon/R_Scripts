Custom_Individual_OptGenMix <- function(max_steps=max_steps, run_removesamples=run_removesamples,
                                        dms=dms, gt_sw_comp=gt_sw_comp, max_t=max_t_indv, N_t_vec=N_t_vec, mvalues=mvalues, ncpu=ncpu, 
                                        unlimited_mvals=unlimited_mvals, measurevals=measurevals, 
                                        samples_to_force=samples_to_force, initial_weights=initial_weights, weights_min=weights_min,
                                        pMAC_mode=pMAC_mode, site_col_name=site_col_name, i_sw_common=i_sw_common, i_sw_rare=i_sw_rare, 
                                        OGM_dir=OGM_dir, threshold_maf=threshold_maf, auto_nt=auto_nt, samples_to_exclude=samples_to_exclude, kinall=NULL){
  
  

  ####Step 1####
  #How many samples do you need to have representative collections, and which individuals should I sample to optimise both rare and common allele capture?
  #How many samples should you optimise for? 
  
  max_steps_random <- max_steps # how many randomisations whould we run? - sometyimes we want this smaller than max_steps for optimisation
  i_ub <- c(1:nrow(gt_sw_comp))
  allvals_common <- mat.or.vec(max_steps_random, length(N_t_vec))
  allvals_rare <- mat.or.vec(max_steps_random, length(N_t_vec))
  rvals <- c()
  for ( i in 1:length(N_t_vec)) {
    iNt <- N_t_vec[i]
    ivals_common <- c()
    ivals_rare <- c()
    cat("\n Running ", i, " ...", iNt, "samples \n")
    for (j in 1:max_steps_random) {
      ran_vec <- rep(0, nrow(gt_sw_comp))
      forcedsamps <- NULL
      excludesamps <- NULL
      if(!is.null(samples_to_force)){
        forcedsamps <- which(rownames(gt_sw_comp)%in%samples_to_force)
        ran_vec <- replace(ran_vec,forcedsamps,1)
      }
      i_ub2 <- which(!i_ub%in%forcedsamps)
      if(!is.null(samples_to_exclude)){
        excludesamps <- which(rownames(gt_sw_comp)%in%samples_to_exclude)
        i_ub2 <- i_ub2[which(!i_ub2%in%excludesamps)]
      }
      #now randomly sample additional samples ontop of what samples are forced and excluded
      ran_vec[sample(i_ub2)[0:(iNt-length(samples_to_force))]] <- 1
      
      common_alleles  <- common_allele_count(gt_sw_comp, ran_vec)
      ivals_common[j] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_common))
      ivals_rare[j] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_rare))
    }
    if (length(N_t_vec)>1){
      allvals_common[,i] <- ivals_common/length(i_sw_common)
      allvals_rare[,i] <- ivals_rare/length(i_sw_rare)
    } else {
      allvals_common <- ivals_common/length(i_sw_common)
      allvals_rare <- ivals_rare/length(i_sw_rare)
    }
  }
  
  allvals_common <- data.frame(allvals_common)
  colnames(allvals_common)[1] <- "X1"
  allvals_common$MAF <- paste0("Common")
  
  allvals_rare <- data.frame(allvals_rare)
  colnames(allvals_rare)[1] <- "X1"
  allvals_rare$MAF <- paste0("Rare")
  
  
  #allvals <- rbind(allvals_common, allvals_common_5pecent,  allvals_common_2pecent, allvals_rare, allvals_rare_5pecent, allvals_rare_2pecent)
  allvals <- rbind(allvals_common, allvals_rare)
  allvals <- data.frame(allvals)
  colnames(allvals) <- c(N_t_vec, "MAF")
  
  allvals2 <- melt(allvals, "MAF")
  allvals2 <- data.frame(allvals2)
  colnames(allvals2) <- c("MAF","nt","prop")
  
  write.table(allvals2, paste0(OGM_dir, "BIGDATA_Randomised_Individuals.csv"), sep=",",quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  
  allvals2min <- data.frame(allvals2 %>% group_by(MAF, nt) %>% slice(which.min(prop)))
  colnames(allvals2min) <- c("MAF","nt","minprop")
  
  
  
  #plot
  minval <- min(c(min(allvals2$prop))) 
  maxval <- max(c(max(allvals2$prop)))
  
  ggplot()+  
    geom_violin(data=allvals2, aes(x = factor(nt), y = prop, group=interaction(nt,MAF), colour=MAF), fill=NA, position = position_dodge(width = 0), scale='width')+ 
    ylim(minval,maxval)+
    theme_minimal()+ 
    geom_hline(yintercept = 0.9, linetype="dashed", alpha=0.5, colour="red")+
    geom_hline(yintercept = 0.5, linetype="dashed", alpha=0.5, colour="blue")+
    #geom_point(data=allvals2, aes(x = factor(nt), y = prop), colour="red", size=0.01) +
    #geom_point(data = Randvals, mapping = aes(x = factor(nt), y = prop),colour=c("black"), size=1)+
    labs(x = "Samples", y = "Allele Proportion", colour="MAF Common Vs Rare")+
    ggtitle(paste0("Total Samples: ",nrow(gt_sw_comp), ". SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position = "right") #c(0.85,0.25))
  # scale_colour_manual(name = "Measure", values = rainbow(length(unique(Optvals$m)))) +   
  # scale_shape_manual(name = "Measure", values = shapeslist[1:length(unique(Optvals$measure))]) 
  
  ggsave(paste0("1. Random Sampling Violin.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  ggplot()+  
    geom_line(data=allvals2min, aes(x = factor(nt), y = minprop, group=MAF, colour=MAF))+ 
    ylim(minval,maxval)+
    theme_minimal()+ 
    geom_hline(yintercept = 0.9, linetype="dashed", alpha=0.5, colour="red")+
    geom_hline(yintercept = 0.5, linetype="dashed", alpha=0.5, colour="blue")+
    #geom_point(data=allvals2, aes(x = factor(nt), y = prop), colour="red", size=0.01) +
    #geom_point(data = Randvals, mapping = aes(x = factor(nt), y = prop),colour=c("black"), size=1)+
    labs(x = "Samples", y = "Allele Proportion", colour="MAF Common Vs Rare")+
    ggtitle(paste0("Maximum AlleleProp Captured from randomisation. Total Samples: ",nrow(gt_sw_comp), ". SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position = "right") #c(0.85,0.25))
  # scale_colour_manual(name = "Measure", values = rainbow(length(unique(Optvals$m)))) +   
  # scale_shape_manual(name = "Measure", values = shapeslist[1:length(unique(Optvals$measure))]) 
  
  ggsave(paste0("2. Random Sampling Mininmmum Line.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  

  
  
  
  #OK ready to optimise based off the randomisation!
  #idenify how many samples to optimsie for in downstream analyses:
  if (auto_nt){
    auto_nt <- allvals2min[which(allvals2min$MAF=="Common" & allvals2min$minprop>0.9),] # find the sample combo where the min random allele prop for 5% MAF reaches over 90% common alleles
    auto_nt <- data.frame(auto_nt[order(auto_nt$nt),])
    auto_nt <- as.numeric(as.character(auto_nt$nt))[1]
    N_t_vec <- auto_nt
  } else {
    N_t_vec <- N_t_vec
  }
  
  max_wts <- rep(1, nrow(gt_sw_comp)) # how many times can an individual be re-sampled? default is only once.
  
  Optvals <- c() # create a new df variable for downstream analyses
  out_alleles <- mat.or.vec(length(N_t_vec),2)
  out_alleles_rare <- mat.or.vec(length(N_t_vec),2)
  allelescapturedfin2 <- c()
  
  for (o in 1:length(measurevals)){
    measure <- measurevals[o]
    ulimM <- unlimited_mvals[o]
    allelescapturedfin <- c()
    set.seed(9825)
    sw_out_list <- list()
    
    for ( i in 1:length(N_t_vec) ) {
      N_t <- N_t_vec[i]
      
      if (mvalues[o] == "auto"){
        m <- N_t #this sets m to the same value as whats being optimised (N_t) rather than the total numebr fo samples in the dataset
      } else {
        m <- mvalues[o]
      }
      
      cat("\n Running ", measure," for ", N_t, "samples ...\n")
      
      #force or exclude any samples?
      if (!is.null(samples_to_force)&&!is.null(samples_to_exclude)){
        forcedsamps <- which(rownames(gt_sw_comp)%in%samples_to_force)
        excludesamps <- which(rownames(gt_sw_comp)%in%samples_to_exclude)
        maxws <- replace(max_wts,forcedsamps,0)
        maxws <- replace(maxws,excludesamps,0)
        max_wts[excludesamps] <- 0
        initial_weights <- propose_initial_weights(nrow(gt_sw_comp), (N_t-length(samples_to_force)), w_max=maxws)
        initial_weights[forcedsamps] <- 1
        weights_min <- rep(0, nrow(gt_sw_comp))
        weights_min <-replace(weights_min,forcedsamps,1)
      } else if (!is.null(samples_to_force)){
        forcedsamps <- which(rownames(gt_sw_comp)%in%samples_to_force)
        maxws <- replace(max_wts,forcedsamps,0)
        initial_weights <- propose_initial_weights(nrow(gt_sw_comp), (N_t-length(samples_to_force)), w_max=maxws)
        initial_weights[forcedsamps] <- 1
        weights_min <- rep(0, nrow(gt_sw_comp))
        weights_min <-replace(weights_min,forcedsamps,1)
      } else if (!is.null(samples_to_exclude)){
        excludesamps <- which(rownames(gt_sw_comp)%in%samples_to_exclude)
        max_wts[excludesamps] <- 0
        initial_weights <- propose_initial_weights(nrow(gt_sw_comp), N_t, w_max=max_wts)
      } 
      
  
      #now run the actual psfs optimisation

      CommOnly <- gt_sw_comp[,which(colnames(gt_sw_comp)%in%rownames(data.frame(i_sw_common)))]
      gt_sw_compComm <- gt_to_minor_alleles(CommOnly)
      
      #RareOnly <- gt_sw_comp[,which(colnames(gt_sw_comp)%in%rownames(data.frame(i_sw_rare_5pecent)))]
      #gt_sw_compRare <- gt_to_minor_alleles(RareOnly)
      
      cat("start optimising common alleles!")
      opt_results <- optimize_single_objective(gt=gt_sw_compComm, sm = NULL, N_t=N_t, measure=measure, max_steps=max_steps, max_t=max_t, m=m, p_depends_delta=FALSE, q=NULL, ncpu=ncpu, weights_max = max_wts,initial_weights = initial_weights, weights_min= weights_min, unlim_m = ulimM)
      
      #cat("start multi optimisation with common and rare alleles!")
      #opt_results <-   optimize_multi_objective(v1=gt_sw_compComm, v2=gt_sw_compRare,  N_t = N_t, measure_1 = measure, measure_2 = measure,  max_steps = max_steps, max_t = max_t,  m = m, p_depends_delta = FALSE, q = NULL, ncpu = ncpu, weights_max = max_wts,initial_weights = initial_weights, weights_min = weights_min, unlim_m = ulimM,  pMAC_mode = pMAC_mode,  kinall=kinall)   
      
      
      
      sw_out_list[[ i ]] <- list(N_t=N_t, m=m, d_opt=opt_results)
      
      OGM_dir_temp <- paste0(OGM_dir,"Temperature_Plots/")
      if (!dir.exists(OGM_dir_temp)) {
        dir.create(OGM_dir_temp, recursive = TRUE)
        cat(paste("Created directory:", OGM_dir_temp, "\n"))
      } else {}
      
      tiff(paste0(OGM_dir_temp ,N_t,"samples ", species, " Temperature Plot T=", max_t, IncludeNA,measure,"m=", m, ".tiff"),
           units = "in", width = 16, height = 10, res = 100)
      par(mfrow = c(1, 1))
      plot(sw_out_list[[i]]$d_opt$value, main= paste0(N_t, "samples T_max = ", max_t," value_1"))  
      dev.off()
      
      #tiff(paste0(OGM_dir_temp ,N_t,"samples ", species, " Temperature Plot T=", max_t, IncludeNA,measure,"m=", m, ".tiff"),
      #     units = "in", width = 16, height = 10, res = 100)
      #  par(mfrow = c(1, 1))
      #  plot(sw_out_list[[i]]$d_opt$value_1, main= paste0(N_t, "samples T_max = ", max_t," value_1"))  
      #  #dev.off()
      
      #tiff(paste0(OGM_dir_temp ,N_t,"samples ", species, " Temperature Plot T=", max_t, IncludeNA,measure,"m=", m, ".tiff"),
      #     units = "in", width = 16, height = 10, res = 100)
      #  par(mfrow = c(1, 1))
      #  plot(sw_out_list[[i]]$d_opt$value_2, main= paste0(N_t, "samples T_max = ", max_t," value_2"))  
      #  dev.off()
      
      
      #remove n samples from the optimized combination to see how allele proportion varies)
      if (run_removesamples==TRUE){
        cat("...running section to remove 1-5 samples from optimised combination")
        for(f in 1:5){
          allelescaptured <- data.frame(matrix(nrow=ncol(combn(N_t_vec[i],f)), ncol=5))
          colnames(allelescaptured) <- c("measure", "m", "nsamps2remove", "nt", "vals")
          for (x in 1:ncol(combn(N_t_vec[i],f))){
            sol_vec2 <- sw_out_list[[i]]$d_opt$weight[max_steps,]
            sol_vec2[which(sw_out_list[[i]]$d_opt$weight[max_steps,]>0)][combn(N_t_vec[i],f)[,x]] <- 0
            common_alleles  <- common_allele_count(gt_sw_comp, sol_vec2)
            allelescaptured$measure[x] <- measure
            allelescaptured$m[x] <- m
            allelescaptured$nsamps2remove[x] <- f
            allelescaptured$nt[x] <- N_t_vec[i]
            allelescaptured$vals_common[x] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_common))/length(i_sw_common)
            allelescaptured$vals_rare[x] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_rare))/length(i_sw_rare)
            cat(paste0(N_t_vec[i],".",x," "))
          }
          allelescapturedfin <- rbind(allelescapturedfin,allelescaptured)
        }
      }
    }
    
    if (run_removesamples==TRUE){
      allelescapturedfin2 <- rbind(allelescapturedfin2,allelescapturedfin)
      samps2removefin <- allelescapturedfin2
    }
    
    
    
    ##produce table
    solution_table <- mat.or.vec(nrow(gt_sw_comp), length(N_t_vec)+4)
    for (i in 1:length(N_t_vec)) { 
      #for each N_t, get the following:
      # 1. get solution table
      solution_table[,i] <- sw_out_list[[i]]$d_opt$weight[max_steps,] 
      
      # 2. find allele prop using common_alleles
      sol_vec <- sw_out_list[[i]]$d_opt$weight[max_steps,]
      
      #plot Super Common Allele frequency for each optimsied sample      
      AllelesInSamps <- gt_sw_comp[which(sol_vec>0),]
      commonALinSampsnum <- which(colnames(AllelesInSamps)%in%rownames(data.frame(i_sw_common)))
      commonALinSamps <- AllelesInSamps[,commonALinSampsnum]
      CAC <- data.frame(common_allele_count(commonALinSamps))
      CAC$propindv <- CAC$minor_allele_counts/length(which(sol_vec>0))
      ggplot()+
        geom_bar(data=CAC, aes(x= minor_allele_counts))+
        theme_minimal()+
        ylab("Frequency")+
        scale_x_continuous("minor_allele_counts", labels = as.character(CAC$minor_allele_counts), breaks = CAC$minor_allele_counts)+
        ggtitle(paste0("Super Common Alleles - common allele count of common alleles \nNumber of common alleles (MAF5%): ", ncol(commonALinSamps), "\n", "# of CA found in less than 10% of select samples: ", length(which(CAC$propindv<0.10)), " (",round(length(which(CAC$propindv<0.10))/length(CAC$propindv)*100, 3), "%)"))
      
      ggsave(paste0("3a. Super Common Alleles in Optimsised Samples.tiff"), path = paste0(OGM_dir), width = 12, height = 8, dpi = 300, units = "in")
      

      common_alleles  <- common_allele_count(gt_sw_comp, sol_vec) #returns: number_common_alleles=number_common_alleles, minor_allele_counts=minor_allele_counts  #common_alleles[[2]]: minor allele count is greater than zero and alleles with a minimum allele freq greater than  0.03 #ie, this asks, which loci were common (> 0.02) in the whole population, and are also represented by two alleles in the proposed conservation population...
      out_alleles[i,1] <- i
      out_alleles[i,2] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_common))
      
      out_alleles_rare[i,1] <- i
      out_alleles_rare[i,2] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_rare))
    }
    
    solution_table[,i+1] <- dms$sample_names
    solution_table[,i+2] <- dms$meta$analyses[,site_col_name]
    solution_table[,i+3] <- dms$meta$lat
    solution_table[,i+4] <- dms$meta$long
    colnames(solution_table) <- c(N_t_vec, "sample", "site", "lat", "long")  #remember to give column names
   
    Optvals <- rbind(Optvals, data.frame(measure= measure,
                                         m=m, 
                                         nt = c(N_t_vec), 
                                         prop_common = cbind(c(out_alleles[,2]/length(i_sw_common))), 
                                         prop_rare = cbind(c(out_alleles_rare[,2]/length(i_sw_rare)))
                                         ))
    
  }
  
  colnames(Optvals) <- c("measure", "m", "nt", "Common", "Rare")
  Optvalsfin <- Optvals 
  Optvals_extend <- matrix(ncol=5, nrow=(nrow(solution_table)-nrow(Optvals)))
  colnames(Optvals_extend) <- colnames(Optvals)
  Optvals_extend <- rbind(Optvals, Optvals_extend)
  Optvals_extend[is.na(Optvals_extend)] <- ""
  
  solution_tablefin <- cbind(solution_table, Optvals_extend)
  
  write.table(solution_tablefin, paste0(OGM_dir, "Optimised_Individual_Combination.csv"), sep=",",quote=FALSE, row.names=FALSE, col.names=TRUE)


  
  #save the range of allele proportion captured when removing samples from optimised combinations
  if (run_removesamples==TRUE){
    rm_sample_min <- data.frame(allelescapturedfin2 %>% group_by(nsamps2remove , nt) %>% slice(which.min(vals_common_5pecent)))
    rm_sample_min$MinMax <- "Min Value"
    rm_sample_max <- data.frame(allelescapturedfin2 %>% group_by(nsamps2remove , nt) %>% slice(which.max(vals_common_5pecent)))
    rm_sample_max$MinMax <- "Max Value"
    rm_sample_range <- rbind(rm_sample_min, rm_sample_max)
    
    rm_sample_range <- rm_sample_range[order(rm_sample_range$MinMax),]
    rm_sample_range <- rm_sample_range[order(rm_sample_range$nt),]
    
    write.csv(rm_sample_range, paste0(OGM_dir, "Range of allele capture removing samples from optimsied combo.csv"),quote=FALSE)
  }
  
  
  #now plot the results!
  library(ggnewscale)
  #optimised vs random
  RandomSamps <- allvals2[which(allvals2$nt==N_t_vec),]
  minval <- min(c(min(RandomSamps$prop), min(Optvalsfin$Rare))) # find he minimum value to constrain plot
  maxval <- max(c(max(RandomSamps$prop), max(Optvalsfin$Common))) # find he minimum value to constrain plot
  shapeslist <- rep(c(15,16,17,18), 10) # get some shapes!
  
  RandomSamps <- RandomSamps[order(RandomSamps$MAF),]
  Optvalsfin <- Optvalsfin[order(Optvalsfin$Common),]
  
  # #Common and Rare - Optimised Vs Random - with all MAF values
  # ggplot()+       
  #   geom_violin(data=RandomSamps, mapping = aes(x = factor(nt), y = prop, colour = MAF), 
  #               fill=NA,size=2, alpha=0.5, position=position_dodge(width=0), scale='width')+
  #   scale_colour_manual(values=c(alpha(rainbow_hcl(6),0.8)))+
  #   new_scale_color()+
  #   geom_point(data = Optvalsfin, mapping = aes(x = factor(nt), y = prop, colour= MAF), 
  #              fill= NA ,size=2, alpha=1, position=position_dodge(width=0))+
  #   scale_colour_manual(values=c(rainbow_hcl(6)))+
  #   labs(x = "MAF", y = "Allele Proportion", colour="MAF", colour="Allele Capture")+
  #   ylim(minval,maxval)+
  #   theme_minimal()+ 
  #   ggtitle(paste0("Common and Rare AlleleProp using 5% MAF VS Random. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
  #   theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  # ggsave(paste0(species, site_col_name,"_Optmised_Vs_Random_Common_and_Rare", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  # 
  #Common Only - Optimised Vs Random - with all MAF values
  # RandomSampsCommonOnly <- RandomSamps[which(RandomSamps$MAF==paste0("1. ", threshold_maf," Common")|RandomSamps$MAF=="2. 5% Common"|RandomSamps$MAF=="3. 2% Common"),]
  # OptvalsfinCommonOnly <- Optvalsfin[which(Optvalsfin$MAF==paste0("1. ", threshold_maf," Common")|Optvalsfin$MAF=="2. 5% Common"|Optvalsfin$MAF=="3. 2% Common"),]
  # minval <- min(c(min(RandomSampsCommonOnly$prop), min(OptvalsfinCommonOnly$prop))) # find he minimum value to constrain plot
  # maxval <- max(c(max(RandomSampsCommonOnly$prop), max(OptvalsfinCommonOnly$prop))) # find he minimum value to constrain plot
  # 
  # ggplot()+       
  #   geom_violin(data=RandomSampsCommonOnly, mapping = aes(x = factor(nt), y = prop, colour = MAF), 
  #               fill=NA,size=2, alpha=0.5, position=position_dodge(width=0), scale="width")+
  #   scale_colour_manual(values=c(alpha(rainbow_hcl(3),0.8)))+
  #   new_scale_color()+
  #   geom_point(data = OptvalsfinCommonOnly, mapping = aes(x = factor(nt), y = prop, colour= MAF), 
  #              fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
  #   scale_colour_manual(values=c(alpha(rainbow_hcl(3),1)))+
  #   labs(x = "N Samples", y = "Allele Proportion", colour="MAF", colour="Allele Capture")+
  #   ylim(minval,maxval)+
  #   theme_minimal()+ 
  #   ggtitle(paste0("Common AlleleProp using 5% MAF VS Random. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
  #   theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  # ggsave(paste0(species, site_col_name,"_Optmised_Vs_Random_Common_Only", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  #Common Only - Optimised Vs Random - with 5% MAF
  RandomSampsCommonOnly <- RandomSamps[which(RandomSamps$MAF=="Common"),]
  minval <- min(c(min(RandomSampsCommonOnly$prop), min(Optvalsfin$Common))) # find he minimum value to constrain plot
  maxval <- max(c(max(RandomSampsCommonOnly$prop), max(Optvalsfin$Common))) # find he minimum value to constrain plot
  
  ggplot()+       
    geom_violin(data=RandomSampsCommonOnly, mapping = aes(x = factor(nt), y = prop, fill=MAF), colour = "black", 
                ,size=2, alpha=0.5, position=position_dodge(width=0), scale="width")+
    geom_point(data = Optvalsfin, mapping = aes(x = factor(nt), y = Common, colour= "Optimised"), 
               fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
    scale_colour_manual(values="red")+
    scale_fill_manual(values=alpha("grey", 0.5))+
    labs(x = "N Samples", y = "Allele Proportion", colour="Optimised", fill="Random")+
    ylim(minval,maxval)+
    theme_minimal()+ 
    ggtitle(paste0("Common AlleleProp using 5% MAF VS Random. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  ggsave(paste0("4. Optmised Vs Random Common Alleles.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
  
  #Rare Only - Optimised Vs Random - with all MAF values
  # RandomSampsRareOnly <- RandomSamps[which(RandomSamps$MAF==paste0("4. ", threshold_maf," Rare")|RandomSamps$MAF=="5. 5% Rare"|RandomSamps$MAF=="6. 2% Rare"),]
  # OptvalsfinRareOnly <- Optvalsfin[which(Optvalsfin$MAF==paste0("4. ", threshold_maf," Rare")|Optvalsfin$MAF=="5. 5% Rare"|Optvalsfin$MAF=="6. 2% Rare"),]
  # minval <- min(c(min(RandomSampsRareOnly$prop), min(OptvalsfinRareOnly$prop))) # find he minimum value to constrain plot
  # maxval <- max(c(max(RandomSampsRareOnly$prop), max(OptvalsfinRareOnly$prop))) # find he minimum value to constrain plot
  # 
  # ggplot()+       
  #   geom_violin(data=RandomSampsRareOnly, mapping = aes(x = factor(nt), y = prop, colour = MAF), 
  #               fill=NA, size=2, alpha=0.5, position=position_dodge(width=0), scale='width')+
  #   scale_colour_manual(values=c(alpha(rainbow_hcl(3),0.8)))+
  #   new_scale_color()+
  #   geom_point(data = OptvalsfinRareOnly, mapping = aes(x = factor(nt), y = prop, colour= MAF), 
  #              fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
  #   scale_colour_manual(values=c(alpha(rainbow_hcl(3),1)))+
  #   labs(x = "N Samples", y = "Allele Proportion", colour="MAF", colour="Allele Capture")+
  #   ylim(minval,maxval)+
  #   theme_minimal()+ 
  #   ggtitle(paste0("Rare AlleleProp using 5% MAF VS Random. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
  #   theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  # ggsave(paste0(species, site_col_name,"_Optmised_Vs_Random_Rare_Only", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
  #Rare Only - Optimised Vs Random - with 5% MAF
  RandomSampsRarenOnly <- RandomSamps[which(RandomSamps$MAF=="Rare"),]
  minval <- min(c(min(RandomSampsRarenOnly$prop), min(Optvalsfin$Rare))) # find he minimum value to constrain plot
  maxval <- max(c(max(RandomSampsRarenOnly$prop), max(Optvalsfin$Rare))) # find he minimum value to constrain plot
  
  ggplot()+       
    geom_violin(data=RandomSampsRarenOnly, mapping = aes(x = factor(nt), y = prop, fill=MAF), colour = "black", 
                ,size=2, alpha=0.5, position=position_dodge(width=0), scale="width")+
    geom_point(data = Optvalsfin, mapping = aes(x = factor(nt), y = Rare, colour= "Optimsed"), 
               fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
    scale_colour_manual(values="red")+
    scale_fill_manual(values=alpha("grey", 0.5))+
    labs(x = "N Samples", y = "Allele Proportion", colour="Optimised", fill="Random")+
    ylim(minval,maxval)+
    theme_minimal()+ 
    ggtitle(paste0("Rare AlleleProp using 5% MAF VS Random. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  ggsave(paste0("5. Optmised Vs Random Rare Alleles.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
  #optimised vs removing samples - Common Only for 5% MAF
  
  if (run_removesamples==TRUE){
    
    # samps2removefinCommonOnly <- samps2removefin[which(samps2removefin$MAF=="vals_common_5pecent"),]
    samps2removefinCommonOnly <- samps2removefin
    minval <- min(c(min(samps2removefinCommonOnly$vals_common_5pecent), min(Optvalsfin$Common))) # find he minimum value to constrain plot
    maxval <- max(c(max(samps2removefinCommonOnly$vals_common_5pecent), max(Optvalsfin$Common))) # find he minimum value to constrain plot
    
    ggplot() +       
      geom_violin(data=samps2removefinCommonOnly, mapping = aes(x = factor(nt), y = vals_common_5pecent, group=interaction(nsamps2remove, nt), colour=factor(nsamps2remove)), 
                  fill=NA,size=2, position=position_dodge(width=0), scale='width')+
      scale_colour_manual(values=c(alpha(rainbow_hcl(length(unique(samps2removefinCommonOnly$nsamps2remove))),0.8)))+
      labs(x = "N Samples", y = "Allele Proportion", colour="n Samps Removed")+
      new_scale_colour()+
      geom_point(data = Optvalsfin, mapping = aes(x = factor(nt), y = Common, group=nt, colour= "Optimised"),  
                 fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
      scale_colour_manual(values=c(alpha("red",1)))+
      labs(x = "N Samples", y = "Allele Proportion", colour="Optimised Combo")+
      ylim(minval,maxval)+
      theme_minimal()+ 
      ggtitle(paste0("Common AlleleProp using 5% MAF removing 1-5 Samples. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
    
    ggsave(paste0("6. Optmised Vs Removing Samples Common Alleles", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
    
    
    #optimised vs removing samples - Rare Only 5% MAF
    # samps2removefinRareOnly <- samps2removefin[which(samps2removefin$MAF=="vals_rare_5pecent"),]
    samps2removefinRareOnly <- samps2removefin
    minval <- min(c(min(samps2removefinRareOnly$vals_rare_5pecent), min(Optvalsfin$Rare))) # find he minimum value to constrain plot
    maxval <- max(c(max(samps2removefinRareOnly$vals_rare_5pecent), max(Optvalsfin$Rare))) # find he minimum value to constrain plot
    
    ggplot() +       
      geom_violin(data=samps2removefinRareOnly, mapping = aes(x = factor(nt), y = vals_rare_5pecent, group=interaction(nsamps2remove, nt), colour=factor(nsamps2remove)), 
                  fill=NA,size=2, position=position_dodge(width=0), scale='width')+
      scale_colour_manual(values=c(alpha(rainbow_hcl(length(unique(samps2removefinRareOnly$nsamps2remove))),0.8)))+
      labs(x = "N Samples", y = "Allele Proportion", colour="n Samps Removed")+
      new_scale_colour()+
      geom_point(data = Optvalsfin, mapping = aes(x = factor(nt), y = Rare, group=nt, colour= "Optimised"), 
                 fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
      scale_colour_manual(values=c(alpha("red",1)))+
      labs(x = "N Samples", y = "Allele Proportion", colour="Optimised Combo",)+
      ylim(minval,maxval)+
      theme_minimal()+ 
      ggtitle(paste0("Rare AlleleProp using 5% MAF removing 1-5 Samples. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
    
    ggsave(paste0("7. Optmised Vs Removing Samples Rare Alleles"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
    
  }
  
  
  
  #### Step 3 ####
  
  #What if the plants at the designated sites aren’t tagged and hard to trace back to one individual?
  #C.	Here’s what you will capture by sampling any tree from the sites and number of trees per site mentioned in the individual optimisation
  
  # read the csv file of the generated optimised solution tables and identify the sites and samples for each
  
  
  solution_table <- data.frame(solution_table)
  allvals_common2 <- mat.or.vec(max_steps_random, length(N_t_vec))
  allvals_rare2 <- mat.or.vec(max_steps_random, length(N_t_vec))

  for (z in 1:length(N_t_vec)) { 
    ran_vec2 <- rep(0, nrow(gt_sw_comp))
    SummaryTab <-c()
    ivals_common2 <- c()
    ivals_rare2 <- c()
    cat("\n Running site/sample variation from optimsied combos ", z, " ...", N_t_vec[z], "samples \n")
    nt_sites <- solution_table$site[which(solution_table[,z]>0)] #identify the sites 
    SummaryTab$site <- unique(nt_sites)
    
    for (g in 1:length(unique(nt_sites))){
      SummaryTab$optsitesnsamps[g] <- length(which(solution_table$site==unique(nt_sites)[g] & solution_table[,z]>0)) #identify how many samples for each site 
    }
    
    for (v in 1:max_steps_random){
      ran_vec2 <- rep(0, nrow(gt_sw_comp))
      
      for (b in 1:length(SummaryTab$site)){
        nswtosamp <- sample(dms$sample_names[which(dms$meta$analyses[,site_col_name]==SummaryTab$site[b])])[1:SummaryTab$optsitesnsamps[b]]
        ran_vec2[which(rownames(gt_sw_comp)%in%c(nswtosamp))] <- 1
      }
      common_alleles  <- common_allele_count(gt_sw_comp, ran_vec2)
      ivals_common2[v] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_common))
      ivals_rare2[v] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_rare))
    }
    
    if (length(N_t_vec)>1){
      allvals_common2[,z] <- ivals_common2/length(i_sw_common)
      allvals_rare2[,z] <- ivals_rare2/length(i_sw_rare)
    } else {
      allvals_common2 <- ivals_common2/length(i_sw_common)
      allvals_rare2 <- ivals_rare2/length(i_sw_rare)
    }  
  }
  
  allvals_common2 <- data.frame(allvals_common2)
  colnames(allvals_common2)[1] <- "X1"
  allvals_common2$MAF <- paste0("Common")
  
  allvals_rare2 <- data.frame(allvals_rare2)
  colnames(allvals_rare2)[1] <- "X1"
  allvals_rare2$MAF <- paste0("Rare")
  
  allvalsver2 <- rbind(allvals_common2, allvals_rare2)
  
  allvalsver2 <- data.frame(allvalsver2)
  colnames(allvalsver2) <- c(N_t_vec, "MAF")
  
  allvalsver22 <- data.frame(melt(allvalsver2, "MAF"))
  colnames(allvalsver22) <- c("MAF","nt","prop")
  
  
  #save the range of allele proportion captured when removing samples from optimised combinations
  rm_sample_min2 <- data.frame(allvalsver22 %>% group_by(MAF , nt) %>% slice(which.min(prop)))
  rm_sample_min2$MinMax <- "Min Value"
  rm_sample_max2 <- data.frame(allvalsver22 %>% group_by(MAF , nt) %>% slice(which.max(prop)))
  rm_sample_max2$MinMax <- "Max Value"
  rm_sample_range2 <- rbind(rm_sample_min2, rm_sample_max2)
  
  rm_sample_range2 <- rm_sample_range2[order(rm_sample_range2$MinMax),]
  rm_sample_range2 <- rm_sample_range2[order(rm_sample_range2$MAF),]
  rm_sample_range2 <- rm_sample_range2[order(rm_sample_range2$nt),]
  
  write.csv(rm_sample_range2, paste0(OGM_dir, "Replacement of optimised samples at the same site.csv"),quote=FALSE)
  
  
  
  #now plot the variation compared to optimised sample combo - with 5% MAF
  
  GeneralSampComboCommonOnly <- allvalsver22[which(allvalsver22$MAF=="Common"),]
  minval <- min(c(min(GeneralSampComboCommonOnly$prop), min(Optvalsfin$Common))) # find he minimum value to constrain plot
  maxval <- max(c(max(GeneralSampComboCommonOnly$prop), max(Optvalsfin$Common))) # find he minimum value to constrain plot
  
  ggplot() +       
    geom_violin(data=GeneralSampComboCommonOnly, mapping = aes(x = factor(nt), y = prop, group=interaction(nt, MAF), colour=factor(nt)), 
                fill=NA,size=2, position=position_dodge(width=0), scale='width')+
    labs(x = "N Samples", y = "Allele Proportion", colour="Samples")+
    scale_colour_manual(values=c(alpha(rainbow_hcl(length(unique(GeneralSampComboCommonOnly$nt))),1)))+
    new_scale_colour()+
    geom_point(data = Optvalsfin, mapping = aes(x = factor(nt), y = Common, group=nt, colour= "Optimised"),  
               fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
    scale_colour_manual(values=c(alpha("red",1)))+
    labs(x = "N Samples", y = "Allele Proportion", colour="Optimised Combo")+
    ylim(minval,maxval)+
    theme_minimal()+ 
    ggtitle(paste0("Common AlleleProp Optimised sampling Vs General sampling across the same sites.5% MAF. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  
  ggsave(paste0("8. Optmised Vs Replacement at the same site Common Alleles.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
  #now plot the rare variation compared to optimised sample combo - with 5% MAF
  
  GeneralSampComboRareOnly <- allvalsver22[which(allvalsver22$MAF=="Rare"),]
  minval <- min(c(min(GeneralSampComboRareOnly$prop), min(Optvalsfin$Rare))) # find he minimum value to constrain plot
  maxval <- max(c(max(GeneralSampComboRareOnly$prop), max(Optvalsfin$Rare))) # find he minimum value to constrain plot
  
  ggplot() +       
    geom_violin(data=GeneralSampComboRareOnly, mapping = aes(x = factor(nt), y = prop, group=interaction(nt, MAF), colour=factor(nt)), 
                fill=NA,size=2, position=position_dodge(width=0), scale='width')+
    labs(x = "N Samples", y = "Allele Proportion", colour="Samples")+
    scale_colour_manual(values=c(alpha(rainbow_hcl(length(unique(GeneralSampComboRareOnly$nt))),1)))+
    new_scale_colour()+
    geom_point(data = Optvalsfin, mapping = aes(x = factor(nt), y = Rare, group=nt, colour= "optimised"),  
               fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
    scale_colour_manual(values=c(alpha("red",1)))+
    labs(x = "N Samples", y = "Allele Proportion", colour="Optimised Combo")+
    ylim(minval,maxval)+
    theme_minimal()+ 
    ggtitle(paste0("Rare AlleleProp Optimised sampling Vs General sampling across the same sites. 5% MAF. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  
  ggsave(paste0("9. Optmised Vs Replacement at the same site Rare Alleles.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
}
