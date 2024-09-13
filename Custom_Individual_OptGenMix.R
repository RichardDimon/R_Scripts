Custom_Individual_OptGenMix <- function(max_steps=max_steps, run_removesamples=run_removesamples,
                                        dms=dms, gt_sw_comp=gt_sw_comp, max_t=max_t, N_t_vec=N_t_vec, mvalues=mvalues, ncpu=ncpu, 
                                        max_wts=max_wts, unlimited_mvals=unlimited_mvals, measurevals=measurevals, 
                                        samples_to_force=samples_to_force, initial_weights=initial_weights, weights_min=weights_min,
                                        pMAC_mode=pMAC_mode, site_col_name=site_col_name, i_sw_common=i_sw_common, i_sw_rare=i_sw_rare, 
                                        i_sw_common_5pecent=i_sw_common_5pecent, i_sw_rare_5pecent=i_sw_rare_5pecent,
                                        i_sw_common_2pecent=i_sw_common_2pecent, i_sw_rare_2pecent=i_sw_rare_2pecent, OGM_dir=OGM_dir,
                                        threshold_maf=threshold_maf){
  
  
  
  ####Step 1####
  #How many samples do you need to have representative collections, and which individuals should I sample to optimise both rare and common allele capture?
  #How many samples should you optimise for? 

  i_ub <- c(1:nrow(gt_sw_comp))
  
  allvals_common <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_rare <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_common_5pecent <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_rare_5pecent <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_common_2pecent <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_rare_2pecent <- mat.or.vec(max_steps, length(N_t_vec))
  rvals <- c()
  for ( i in 1:length(N_t_vec)) {
    iNt <- N_t_vec[i]
    ivals_common <- c()
    ivals_rare <- c()
    ivals_common_5pecent <- c()
    ivals_rare_5pecent <- c()
    ivals_common_2pecent <- c()
    ivals_rare_2pecent <- c()
    cat("\n Running ", i, " ...", iNt, "samples \n")
    for (j in 1:max_steps) {
      ran_vec <- rep(0, nrow(gt_sw_comp))
      ran_vec[sample(i_ub)[1:iNt]] <- 1
      common_alleles  <- common_allele_count(gt_sw_comp, ran_vec)
      ivals_common[j] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_common))
      ivals_rare[j] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_rare))
      ivals_common_5pecent[j] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_common_5pecent))
      ivals_rare_5pecent[j] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_rare_5pecent))
      ivals_common_2pecent[j] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_common_2pecent))
      ivals_rare_2pecent[j] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_rare_2pecent))
    }
    allvals_common[,i] <- ivals_common/length(i_sw_common)
    allvals_rare[,i] <- ivals_rare/length(i_sw_rare)
    allvals_common_5pecent[,i] <- ivals_common_5pecent/length(i_sw_common_5pecent)
    allvals_rare_5pecent[,i] <- ivals_rare_5pecent/length(i_sw_rare_5pecent)
    allvals_common_2pecent[,i] <- ivals_common_2pecent/length(i_sw_common_2pecent)
    allvals_rare_2pecent[,i] <- ivals_rare_2pecent/length(i_sw_rare_2pecent)
  }
  
  allvals_common <- data.frame(allvals_common)
  allvals_rare <- data.frame(allvals_rare)
  allvals_common$MAF <- paste0("1. ", threshold_maf," Common")
  allvals_rare$MAF <- paste0("4. ", threshold_maf," Rare")
  allvals_common_5pecent <- data.frame(allvals_common_5pecent)
  allvals_rare_5pecent <- data.frame(allvals_rare_5pecent)
  allvals_common_5pecent$MAF <- "2. 5% Common"
  allvals_rare_5pecent$MAF <- "5. 5% Rare"
  allvals_common_2pecent <- data.frame(allvals_common_2pecent)
  allvals_rare_2pecent <- data.frame(allvals_rare_2pecent)
  allvals_common_2pecent$MAF <- "3. 2% Common"
  allvals_rare_2pecent$MAF <- "6. 2% Rare"
  allvals <- rbind(allvals_common, allvals_common_5pecent,  allvals_common_2pecent, allvals_rare, allvals_rare_5pecent, allvals_rare_2pecent)
  
  allvals <- data.frame(allvals)
  colnames(allvals) <- c(N_t_vec, "MAF")
  
  allvals2 <- melt(allvals, "MAF")
  allvals2 <- data.frame(allvals2)
  colnames(allvals2) <- c("MAF","nt","prop")
  
  allvals2max <- data.frame(allvals2 %>% group_by(MAF, nt) %>% slice(which.max(prop)))
  colnames(allvals2max) <- c("MAF","nt","minprop")
  
  
  #idenify how many samples to optimsie for in downstream analyses:
  auto_nt <- allvals2max[which(allvals2max$MAF=="3. 2% Common" & allvals2max$minprop>0.9),] # find the sample combo where the max random allele prop for 2&MAF reaches over 90% common alleles
  auto_nt <- data.frame(auto_nt[order(auto_nt$nt),])
  auto_nt <- as.numeric(as.character(auto_nt$nt))[1]
  
  
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
  
  ggsave(paste0("1. ",species, site_col_name,"_Individual_Randomisation_Violin.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  ggplot()+  
    geom_line(data=allvals2max, aes(x = factor(nt), y = minprop, group=MAF, colour=MAF))+ 
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
  
  ggsave(paste0("2. ",species, site_col_name,"_Individual_Maximum_AlleleProp_Capture_from_Randomisation_Line.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
  
  
  
  
  
  
  
  
  
  
  #OK ready to optimise based off the randomisation!
  
  
  N_t_vec <- c(auto_nt-4, auto_nt-2, auto_nt, auto_nt+2, auto_nt+4) #we know this number from previously
  
  
  Optvals <- c() # create a new df variable for downstream analyses

  
  out_alleles <- mat.or.vec(length(N_t_vec),2)
  out_alleles_rare <- mat.or.vec(length(N_t_vec),2)
  out_alleles_5pecent <- mat.or.vec(length(N_t_vec),2)
  out_alleles_rare_5pecent <- mat.or.vec(length(N_t_vec),2)
  out_alleles_2pecent <- mat.or.vec(length(N_t_vec),2)
  out_alleles_rare_2pecent <- mat.or.vec(length(N_t_vec),2)
  
  
  allelescapturedfin2 <- c()
  for (o in 1:length(measurevals)){
    measure <- measurevals[o]
    ulimM <- unlimited_mvals[o]
    m <- mvalues[o]
    allelescapturedfin <- c()
    set.seed(9825)
    sw_out_list <- list()
    tiff(paste0(OGM_dir,"3. ", species, "Temperature Plots T=", max_t, IncludeNA,measure,"m=", m, ".tiff"),
         units = "in", width = 10, height = 14, res = 100)
    par(mfrow = c(length(N_t_vec), 1))
    for ( i in 1:length(N_t_vec) ) {
      N_t <- N_t_vec[i]
      cat("\n Running ", measure," for ", N_t, "samples ...\n")
      
      if (any(samples_to_force)){
        maxws <- replace(max_wts,samples_to_force,0)
        initial_weights <- propose_initial_weights(nrow(gt_sw_comp), (N_t-length(samples_to_force)), w_max=maxws)
        initial_weights[samples_to_force] <- 1
        weights_min <- rep(0, nrow(gt_sw_comp))
        weights_min <-replace(weights_min,samples_to_force,1)
      }
      
      if (measure=="psfs"){ 
        gt_sw_comp2 <- gt_to_minor_alleles(gt_sw_comp) 
        opt_results <- optimize_single_objective(gt=gt_sw_comp2, N_t=N_t, measure=measure, max_steps=max_steps, max_t=max_t, m=m, p_depends_delta=FALSE, q=NULL, ncpu=ncpu, weights_max = max_wts,initial_weights = initial_weights, weights_min= weights_min, unlim_m = ulimM)
      } else {print("uh oh! add more options here from other versions of OptGenMix scripts")}
      
      sw_out_list[[ i ]] <- list(N_t=N_t, m=m, d_opt=opt_results)
      plot(sw_out_list[[i]]$d_opt$value, main= paste0(site_col_name, " T_max = ", max_t," ",N_t, "samps ", IncludeNA, measure)) 
      
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
            allelescaptured$vals_common_5pecent[x] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_common_5pecent))/length(i_sw_common_5pecent)
            allelescaptured$vals_rare_5pecent[x] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_rare_5pecent))/length(i_sw_rare_5pecent)
            allelescaptured$vals_common_2pecent[x] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_common_2pecent))/length(i_sw_common_2pecent)
            allelescaptured$vals_rare_2pecent[x] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_rare_2pecent))/length(i_sw_rare_2pecent)
            cat(paste0(N_t_vec[i],".",x," "))
          }
          allelescapturedfin <- rbind(allelescapturedfin,allelescaptured)
        }
      }
    }
    
    dev.off()
    
    if (run_removesamples==TRUE){
      allelescapturedfin2 <- rbind(allelescapturedfin2,allelescapturedfin)
      samps2removefin <- allelescapturedfin2
    }
    
    ##produce table
    solution_table <- mat.or.vec(nrow(gt_sw_comp2), length(N_t_vec)+4)
    for (i in 1:length(N_t_vec)) { 
      #for each N_t, get the following:
      # 1. get solution table
      solution_table[,i] <- sw_out_list[[i]]$d_opt$weight[max_steps,] 
      
      # 2. find allele prop using common_alleles
      sol_vec <- sw_out_list[[i]]$d_opt$weight[max_steps,]
      common_alleles  <- common_allele_count(gt_sw_comp, sol_vec) #returns: number_common_alleles=number_common_alleles, minor_allele_counts=minor_allele_counts  #common_alleles[[2]]: minor allele count is greater than zero and alleles with a minimum allele freq greater than  0.03 #ie, this asks, which loci were common (> 0.02) in the whole population, and are also represented by two alleles in the proposed conservation population...
      out_alleles[i,1] <- i
      out_alleles[i,2] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_common))
      
      out_alleles_5pecent[i,1] <- i
      out_alleles_5pecent[i,2] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_common_5pecent))
      
      out_alleles_2pecent[i,1] <- i
      out_alleles_2pecent[i,2] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_common_2pecent))
      
      out_alleles_rare[i,1] <- i
      out_alleles_rare[i,2] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_rare))
      
      out_alleles_rare_5pecent[i,1] <- i
      out_alleles_rare_5pecent[i,2] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_rare_5pecent))
      
      out_alleles_rare_2pecent[i,1] <- i
      out_alleles_rare_2pecent[i,2] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_rare_2pecent))
    }
    
    solution_table[,i+1] <- dms$sample_names
    solution_table[,i+2] <- dms$meta$analyses[,site_col_name]
    solution_table[,i+3] <- dms$meta$lat
    solution_table[,i+4] <- dms$meta$long
    colnames(solution_table) <- c(N_t_vec, "sample", "site", "lat", "long")  #remember to give column names
    write.table(solution_table, paste0(OGM_dir, species,"_",site_col_name,"_Individual_Combinations_singleton=",threshold_maf, "_maxt=", max_t,"_",IncludeNA,measure,"m=",m,".csv"), sep=",",quote=FALSE, row.names=FALSE, col.names=TRUE)
    Optvals <- rbind(Optvals, data.frame(measure= measure,
                                         m=m, 
                                         nt = c(N_t_vec), 
                                         prop_common = cbind(c(out_alleles[,2]/length(i_sw_common))), 
                                         prop_rare = cbind(c(out_alleles_rare[,2]/length(i_sw_rare))),
                                         prop_common_5pecent = cbind(c(out_alleles_5pecent[,2]/length(i_sw_common_5pecent))), 
                                         prop_rare_5pecent = cbind(c(out_alleles_rare_5pecent[,2]/length(i_sw_rare_5pecent))),
                                         prop_common_2pecent = cbind(c(out_alleles_2pecent[,2]/length(i_sw_common_2pecent))), 
                                         prop_rare_2pecent = cbind(c(out_alleles_rare_2pecent[,2]/length(i_sw_rare_2pecent)))
                                         
    ))
    
  }
  
  colnames(Optvals) <- c("measure", "m", "nt",  paste0("1. ", threshold_maf," Common"), paste0("4. ", threshold_maf," Rare"), "2. 5% Common", "5. 5% Rare", "3. 2% Common", "6. 2% Rare")
  
  Optvals_AC <- melt(Optvals[3:9], "nt")
  
  Optvalsfin <- Optvals_AC
  # Optvalsfin$nt <- Optvals$nt
  colnames(Optvalsfin) <- c("nt","MAF","prop")
  
  #save optimised allele proportion cpatured for each nt
  write.csv(Optvalsfin, paste0(OGM_dir, species,"_",site_col_name,"_Optimsied common and rare proportion captured across different MAF thresholds.csv"),quote=FALSE)
  
  
  #save the range of allele proportion captured when removing samples from optimised combinations
  if (run_removesamples==TRUE){
    rm_sample_min <- data.frame(allelescapturedfin2 %>% group_by(nsamps2remove , nt) %>% slice(which.min(vals_common_5pecent)))
    rm_sample_min$MinMax <- "Min Value"
    rm_sample_max <- data.frame(allelescapturedfin2 %>% group_by(nsamps2remove , nt) %>% slice(which.max(vals_common_5pecent)))
    rm_sample_max$MinMax <- "Max Value"
    rm_sample_range <- rbind(rm_sample_min, rm_sample_max)
    
    rm_sample_range <- rm_sample_range[order(rm_sample_range$MinMax),]
    rm_sample_range <- rm_sample_range[order(rm_sample_range$nt),]
    
    write.csv(rm_sample_range, paste0(OGM_dir, species,"_",site_col_name,"_range of AlleleProp captured removing samples from optimsied combos.csv"),quote=FALSE)
  }
  
  
  #now plot the results!
  library(ggnewscale)
  #optimised vs random
  RandomSamps <- allvals2[which(allvals2$nt==N_t_vec),]
  minval <- min(c(min(RandomSamps$prop), min(Optvalsfin$prop))) # find he minimum value to constrain plot
  maxval <- max(c(max(RandomSamps$prop), max(Optvalsfin$prop))) # find he minimum value to constrain plot
  shapeslist <- rep(c(15,16,17,18), 10) # get some shapes!
  
  RandomSamps <- RandomSamps[order(RandomSamps$MAF),]
  Optvalsfin$MAF <- as.character(Optvalsfin$MAF)
  Optvalsfin <- Optvalsfin[order(Optvalsfin$MAF),]
  
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
  RandomSampsCommonOnly <- RandomSamps[which(RandomSamps$MAF=="2. 5% Common"),]
  OptvalsfinCommonOnly <- Optvalsfin[which(Optvalsfin$MAF=="2. 5% Common"),]
  minval <- min(c(min(RandomSampsCommonOnly$prop), min(OptvalsfinCommonOnly$prop))) # find he minimum value to constrain plot
  maxval <- max(c(max(RandomSampsCommonOnly$prop), max(OptvalsfinCommonOnly$prop))) # find he minimum value to constrain plot
  
  ggplot()+       
    geom_violin(data=RandomSampsCommonOnly, mapping = aes(x = factor(nt), y = prop, fill=MAF), colour = "black", 
                ,size=2, alpha=0.5, position=position_dodge(width=0), scale="width")+
    geom_point(data = OptvalsfinCommonOnly, mapping = aes(x = factor(nt), y = prop, colour= MAF), 
               fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
    scale_colour_manual(values="red")+
    scale_fill_manual(values=alpha("grey", 0.5))+
    labs(x = "N Samples", y = "Allele Proportion", colour="Optimised", fill="Random")+
    ylim(minval,maxval)+
    theme_minimal()+ 
    ggtitle(paste0("Common AlleleProp using 5% MAF VS Random. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  ggsave(paste0("4. ", species, site_col_name,"_Optmised_Vs_Random_Common_Only_5%_MAF_Only", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
  
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
  RandomSampsRarenOnly <- RandomSamps[which(RandomSamps$MAF=="5. 5% Rare"),]
  OptvalsfinRareOnly <- Optvalsfin[which(Optvalsfin$MAF=="5. 5% Rare"),]
  minval <- min(c(min(RandomSampsRarenOnly$prop), min(OptvalsfinRareOnly$prop))) # find he minimum value to constrain plot
  maxval <- max(c(max(RandomSampsRarenOnly$prop), max(OptvalsfinRareOnly$prop))) # find he minimum value to constrain plot
  
  ggplot()+       
    geom_violin(data=RandomSampsRarenOnly, mapping = aes(x = factor(nt), y = prop, fill=MAF), colour = "black", 
                ,size=2, alpha=0.5, position=position_dodge(width=0), scale="width")+
    geom_point(data = OptvalsfinRareOnly, mapping = aes(x = factor(nt), y = prop, colour= MAF), 
               fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
    scale_colour_manual(values="red")+
    scale_fill_manual(values=alpha("grey", 0.5))+
    labs(x = "N Samples", y = "Allele Proportion", colour="Optimised", fill="Random")+
    ylim(minval,maxval)+
    theme_minimal()+ 
    ggtitle(paste0("Rare AlleleProp using 5% MAF VS Random. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  ggsave(paste0("5. ", species, site_col_name,"_Optmised_Vs_Random_Rare_Only_5%_MAF_Only", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
  #optimised vs removing samples - Common Only for 5% MAF
  
  if (run_removesamples==TRUE){
  
    # samps2removefinCommonOnly <- samps2removefin[which(samps2removefin$MAF=="vals_common_5pecent"),]
    samps2removefinCommonOnly <- samps2removefin
    OptvalsfinCommonOnly <- Optvalsfin[which(Optvalsfin$MAF=="2. 5% Common"),]
    
    minval <- min(c(min(samps2removefinCommonOnly$vals_common_5pecent), min(OptvalsfinCommonOnly$prop))) # find he minimum value to constrain plot
    maxval <- max(c(max(samps2removefinCommonOnly$vals_common_5pecent), max(OptvalsfinCommonOnly$prop))) # find he minimum value to constrain plot
    
    ggplot() +       
      geom_violin(data=samps2removefinCommonOnly, mapping = aes(x = factor(nt), y = vals_common_5pecent, group=interaction(nsamps2remove, nt), colour=factor(nsamps2remove)), 
                  fill=NA,size=2, position=position_dodge(width=0), scale='width')+
      scale_colour_manual(values=c(alpha(rainbow_hcl(length(unique(samps2removefinCommonOnly$nsamps2remove))),0.8)))+
      labs(x = "N Samples", y = "Allele Proportion", colour="n Samps Removed")+
      new_scale_colour()+
      geom_point(data = OptvalsfinCommonOnly, mapping = aes(x = factor(nt), y = prop, group=interaction(nt, MAF), colour= MAF),  
                 fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
      scale_colour_manual(values=c(alpha("red",1)))+
      labs(x = "N Samples", y = "Allele Proportion", colour="Optimised Combo")+
      ylim(minval,maxval)+
      theme_minimal()+ 
      ggtitle(paste0("Common AlleleProp using 5% MAF removing 1-5 Samples. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
    
    ggsave(paste0("6. ", species, site_col_name,"_Optmised_Vs_RemovingSamples_Common_Only", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
    
    
    #optimised vs removing samples - Rare Only 5% MAF
    # samps2removefinRareOnly <- samps2removefin[which(samps2removefin$MAF=="vals_rare_5pecent"),]
    samps2removefinRareOnly <- samps2removefin
    OptvalsfinRareOnly <- Optvalsfin[which(Optvalsfin$MAF=="5. 5% Rare"),]
    minval <- min(c(min(samps2removefinRareOnly$vals_rare_5pecent), min(OptvalsfinRareOnly$prop))) # find he minimum value to constrain plot
    maxval <- max(c(max(samps2removefinRareOnly$vals_rare_5pecent), max(OptvalsfinRareOnly$prop))) # find he minimum value to constrain plot
    
    ggplot() +       
      geom_violin(data=samps2removefinRareOnly, mapping = aes(x = factor(nt), y = vals_rare_5pecent, group=interaction(nsamps2remove, nt), colour=factor(nsamps2remove)), 
                  fill=NA,size=2, position=position_dodge(width=0), scale='width')+
      scale_colour_manual(values=c(alpha(rainbow_hcl(length(unique(samps2removefinRareOnly$nsamps2remove))),0.8)))+
      labs(x = "N Samples", y = "Allele Proportion", colour="n Samps Removed")+
      new_scale_colour()+
      geom_point(data = OptvalsfinRareOnly, mapping = aes(x = factor(nt), y = prop, group=interaction(nt, MAF), colour= MAF), 
                 fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
      scale_colour_manual(values=c(alpha("red",1)))+
      labs(x = "N Samples", y = "Allele Proportion", colour="Optimised Combo",)+
      ylim(minval,maxval)+
      theme_minimal()+ 
      ggtitle(paste0("Rare AlleleProp using 5% MAF removing 1-5 Samples. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
    
    ggsave(paste0("7. ",species, site_col_name,"_Optmised_Vs_RemovingSamples_Rare_Only", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
    
  }
  
  
  
  #### Step 3 ####
  
  #What if the plants at the designated sites aren’t tagged and hard to trace back to one individual?
  #C.	Here’s what you will capture by sampling any tree from the sites and number of trees per site mentioned in the individual optimisation
  
  # read the csv file of the generated optimised solution tables and identify the sites and samples for each
  
  
  solution_table <- data.frame(solution_table)
  allvals_common2 <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_rare2 <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_common_5pecent2 <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_rare_5pecent2 <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_common_2pecent2 <- mat.or.vec(max_steps, length(N_t_vec))
  allvals_rare_2pecent2 <- mat.or.vec(max_steps, length(N_t_vec))
  
  
  
  for (z in 1:length(N_t_vec)) { 
    ran_vec2 <- rep(0, nrow(gt_sw_comp))
    SummaryTab <-c()
    ivals_common2 <- c()
    ivals_rare2 <- c()
    ivals_common_5pecent2 <- c()
    ivals_rare_5pecent2 <- c()
    ivals_common_2pecent2 <- c()
    ivals_rare_2pecent2 <- c()
    cat("\n Running site/sample variation from optimsied combos ", z, " ...", N_t_vec[z], "samples \n")
    nt_sites <- solution_table$site[which(solution_table[,z]>0)] #identify the sites 
    SummaryTab$site <- unique(nt_sites)
    
    for (g in 1:length(unique(nt_sites))){
      SummaryTab$optsitesnsamps[g] <- length(which(solution_table$site==unique(nt_sites)[g] & solution_table[,z]>0)) #identify how many samples for each site 
    }
    
    for (v in 1:max_steps){
      ran_vec2 <- rep(0, nrow(gt_sw_comp))
      
      for (b in 1:length(SummaryTab$site)){
        nswtosamp <- sample(dms$sample_names[which(dms$meta$analyses[,site_col_name]==SummaryTab$site[b])])[1:SummaryTab$optsitesnsamps[b]]
        ran_vec2[which(rownames(gt_sw_comp)%in%c(nswtosamp))] <- 1
      }
      common_alleles  <- common_allele_count(gt_sw_comp, ran_vec2)
      ivals_common2[v] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_common))
      ivals_rare2[v] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_rare))
      ivals_common_5pecent2[v] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_common_5pecent))
      ivals_rare_5pecent2[v] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_rare_5pecent))
      ivals_common_2pecent2[v] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_common_2pecent))
      ivals_rare_2pecent2[v] <- length( intersect( which(common_alleles[[2]] > 0), i_sw_rare_2pecent))
    }
    
    allvals_common2[,z] <- ivals_common2/length(i_sw_common)
    allvals_rare2[,z] <- ivals_rare2/length(i_sw_rare)
    allvals_common_5pecent2[,z] <- ivals_common_5pecent2/length(i_sw_common_5pecent)
    allvals_rare_5pecent2[,z] <- ivals_rare_5pecent2/length(i_sw_rare_5pecent)
    allvals_common_2pecent2[,z] <- ivals_common_2pecent2/length(i_sw_common_2pecent)
    allvals_rare_2pecent2[,z] <- ivals_rare_2pecent2/length(i_sw_rare_2pecent)
  }
  
  
  allvals_common2 <- data.frame(allvals_common2)
  allvals_rare2 <- data.frame(allvals_rare2)
  allvals_common2$MAF <- paste0("1. ", threshold_maf," Common")
  allvals_rare2$MAF <- paste0("4. ", threshold_maf," Rare")
  allvals_common_5pecent2 <- data.frame(allvals_common_5pecent2)
  allvals_rare_5pecent2 <- data.frame(allvals_rare_5pecent2)
  allvals_common_5pecent2$MAF <- "2. 5% Common"
  allvals_rare_5pecent2$MAF <- "5. 5% Rare"
  allvals_common_2pecent2 <- data.frame(allvals_common_2pecent2)
  allvals_rare_2pecent2 <- data.frame(allvals_rare_2pecent2)
  allvals_common_2pecent2$MAF <- "3. 2% Common"
  allvals_rare_2pecent2$MAF <- "6. 2% Rare"
  allvalsver2 <- rbind(allvals_common2, allvals_common_5pecent2,  allvals_common_2pecent2, allvals_rare2, allvals_rare_5pecent2, allvals_rare_2pecent2)
  
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
  
  write.csv(rm_sample_range2, paste0(OGM_dir, species,"_",site_col_name,"_range of AlleleProp captured when considering sampling optimised combos from the general sites.csv"),quote=FALSE)
  
  
  
  #now plot the variation compared to optimised sample combo - with 5% MAF
  
  GeneralSampComboCommonOnly <- allvalsver22[which(allvalsver22$MAF=="2. 5% Common"),]
  OptvalsfinCommonOnly <- Optvalsfin[which(Optvalsfin$MAF=="2. 5% Common"),]
  
  minval <- min(c(min(GeneralSampComboCommonOnly$prop), min(OptvalsfinCommonOnly$prop))) # find he minimum value to constrain plot
  maxval <- max(c(max(GeneralSampComboCommonOnly$prop), max(OptvalsfinCommonOnly$prop))) # find he minimum value to constrain plot
  
  ggplot() +       
    geom_violin(data=GeneralSampComboCommonOnly, mapping = aes(x = factor(nt), y = prop, group=interaction(nt, MAF), colour=factor(nt)), 
                fill=NA,size=2, position=position_dodge(width=0), scale='width')+
    labs(x = "N Samples", y = "Allele Proportion", colour="Samples")+
    scale_colour_manual(values=c(alpha(rainbow_hcl(length(unique(GeneralSampComboCommonOnly$nt))),1)))+
    new_scale_colour()+
    geom_point(data = OptvalsfinCommonOnly, mapping = aes(x = factor(nt), y = prop, group=interaction(nt, MAF), colour= MAF),  
               fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
    scale_colour_manual(values=c(alpha("red",1)))+
    labs(x = "N Samples", y = "Allele Proportion", colour="Optimised Combo")+
    ylim(minval,maxval)+
    theme_minimal()+ 
    ggtitle(paste0("Common AlleleProp Optimised sampling Vs General sampling across the same sites.5% MAF. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  
  ggsave(paste0("8. ", species, site_col_name,"_Optmised_Vs_General_sampling_across_same_sites_Common_Only", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
  #now plot the rare variation compared to optimised sample combo - with 5% MAF
  
  GeneralSampComboRareOnly <- allvalsver22[which(allvalsver22$MAF=="5. 5% Rare"),]
  OptvalsfinRareOnly <- Optvalsfin[which(Optvalsfin$MAF=="5. 5% Rare"),]
  
  minval <- min(c(min(GeneralSampComboRareOnly$prop), min(OptvalsfinRareOnly$prop))) # find he minimum value to constrain plot
  maxval <- max(c(max(GeneralSampComboRareOnly$prop), max(OptvalsfinRareOnly$prop))) # find he minimum value to constrain plot
  
  ggplot() +       
    geom_violin(data=GeneralSampComboRareOnly, mapping = aes(x = factor(nt), y = prop, group=interaction(nt, MAF), colour=factor(nt)), 
                fill=NA,size=2, position=position_dodge(width=0), scale='width')+
    labs(x = "N Samples", y = "Allele Proportion", colour="Samples")+
    scale_colour_manual(values=c(alpha(rainbow_hcl(length(unique(GeneralSampComboRareOnly$nt))),1)))+
    new_scale_colour()+
    geom_point(data = OptvalsfinRareOnly, mapping = aes(x = factor(nt), y = prop, group=interaction(nt, MAF), colour= MAF),  
               fill=NA,size=2, alpha=1, position=position_dodge(width=0))+
    scale_colour_manual(values=c(alpha("red",1)))+
    labs(x = "N Samples", y = "Allele Proportion", colour="Optimised Combo")+
    ylim(minval,maxval)+
    theme_minimal()+ 
    ggtitle(paste0("Rare AlleleProp Optimised sampling Vs General sampling across the same sites. 5% MAF. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
    theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
  
  ggsave(paste0("9. ", species, site_col_name,"_Optmised_Vs_General_sampling_across_same_sites_Rare_Only", max_t,"_",IncludeNA,".tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
  
  
}
  
  