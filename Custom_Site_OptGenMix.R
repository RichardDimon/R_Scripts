Custom_Site_OptGenMix <- function(max_steps=max_steps,samplethreshold=samplethreshold, sitethreshold=sitethreshold,
                                  dms=dms, gt_sw_comp=gt_sw_comp, max_t=max_t_sites, mvalues=mvalues, ncpu=ncpu, 
                                  unlimited_mvals=unlimited_mvals, measurevals=measurevals, 
                                  sites_to_force=sites_to_force, initial_weights=initial_weights, weights_min=weights_min,
                                  pMAC_mode=pMAC_mode, site_col_name=site_col_name, i_sw_common=i_sw_common, i_sw_rare=i_sw_rare, 
                                  OGM_dir=OGM_dir, threshold_maf=threshold_maf, manual_sites=manual_sites, manual_sampspersite=manual_sampspersite,
                                  sites_to_exclude=sites_to_exclude, kinall=NULL){
  
  library(ggplot2)
  library(ggnewscale)
  #remove any sites that have less than 5 samples
  sites_high_missing_samples_removed <- table(dms$meta$analyses[,site_col_name])
  not_n5_sites <- as.vector(names(which(sites_high_missing_samples_removed<samplethreshold)))
  not_n5_samples <- dms$sample_names[which(!(dms$meta$analyses[,site_col_name] %in% not_n5_sites))]
  dms <- remove.by.list(dms, not_n5_samples)
  
  cat("\n Removing the following sites as they don't meet the provided sample threshold of ", samplethreshold, " samples: \n")
  print(table(dms$meta$analyses[,site_col_name][which(dms$meta$analyses[,site_col_name] %in% not_n5_sites)]))
  
  cat("\n The remaining sites for site optimsiation are: \n")
  print(table(dms$meta$analyses[,site_col_name]))
  
  if (any(sites_to_force%in%not_n5_sites)){
    cat("\n Awwwww, SNAP! some of the sites provided in the forced site list have been removed! change your forced site list and try again :) ")
  } else if (any(sites_to_exclude%in%not_n5_sites)){
    cat("NOTE: some of the sites you have provided in your site exclude list have been removed! consider changing this list")
  } else {
    if (auto_nt) {
      ##### Common allele proportion for entire Genetic neighborhood for 5 samples per site #####
      
      RandomSiteIterations <- max_steps #may not need10,000 loops for each site combo as it takes too long, 1000 may be plenty for most analyses. if so change this here
      cat("calculating common and rare alleles across differnt MAF theshold for", samplethreshold, " samples, and ", sitethreshold, " sites for ", RandomSiteIterations, " iterations")
    
      All5Sites_MAF5 <- Common_Allele_Prop_Random_Sites(dms = dms, gt_sw_comp=gt_sw_comp, analysis=site_col_name, NumSteps=RandomSiteIterations, i_sw_common=i_sw_common, samplethreshold = samplethreshold, sitethreshold = sitethreshold, sites_to_exclude=sites_to_exclude, sites_to_force=sites_to_force)
      All5Sites_rare_MAF5 <- Common_Allele_Prop_Random_Sites(dms = dms, gt_sw_comp=gt_sw_comp, analysis=site_col_name, NumSteps=RandomSiteIterations, i_sw_common=i_sw_rare, samplethreshold = samplethreshold, sitethreshold = sitethreshold, sites_to_exclude=sites_to_exclude, sites_to_force=sites_to_force)
      
      All5Sites_MAF5$Group <- "5% MAF"
      All5Sites_MAF5$Group2 <- "Common"
      findata <- All5Sites_MAF5
      

      All5Sites_rare_MAF5$Group <- "5% MAF"
      All5Sites_rare_MAF5$Group2 <- "rare"
      findata_rare <- All5Sites_rare_MAF5
      
      findata_common_and_rare <- rbind(findata, findata_rare)
      write.table(findata_common_and_rare, paste0(OGM_dir, "1. BIGDATA Randomised Sites.csv"), sep=",",quote=FALSE, row.names=FALSE, col.names=TRUE)
      
      allvals2minsite <- data.frame(findata_common_and_rare %>% group_by(t_num_indv , n_sites_sel, Group, Group2) %>% slice(which.min(Allele )))
      
      
      
      #plot Common and Rare
      ggplot() +
        geom_violin(data=findata_common_and_rare, aes(x=interaction(n_sites_sel,t_num_indv),y=Allele, colour=interaction(Group,Group2)),
                    fill=NA, position = position_dodge(width = 0), scale="width")+
        
        scale_x_discrete(labels= unique(paste0(findata_common_and_rare$t_num_indv," individuals\n(", findata_common_and_rare$indv_p_site," from ",findata_common_and_rare$n_sites_sel," sites)")))+
        theme_bw()+
        labs(fill= "Common Vs Rare", title=paste0(species," sampling strategy"), 
             y="Alleles (proportion)", x="Number of Sites/ \n Total Samples", colour="")+
        geom_hline(yintercept = 0.9, linetype="dashed", alpha=0.5, colour="red")+
        geom_hline(yintercept = 0.5, linetype="dashed", alpha=0.5, colour="blue")+
        #facet_nested(~t_num_indv+n_sites_sel, scales = "free_x",  space = "free", switch="x")+theme_classic()+
        theme(panel.spacing=unit(0,"lines"), axis.title.x = element_blank(),
              axis.text = element_text(size=8, face="bold"),
              strip.background=element_rect(color="grey30", fill="grey90"), 
              panel.grid = element_blank(), legend.position = "bottom",
              legend.title = element_blank(),
              legend.background = element_rect(fill=alpha("grey", 0.2)))
      
      ggsave(paste0("1. Site Randomisation Violin.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
      
      # ggplot() +
      #   geom_line(data=allvals2minsite, aes(x=interaction(n_sites_sel,t_num_indv),y=Allele, group=interaction(Group, Group2), colour=interaction(Group, Group2)))+
      #   scale_x_discrete(labels= unique(paste0(findata_common_and_rare$t_num_indv," individuals\n(", findata_common_and_rare$indv_p_site," from ",findata_common_and_rare$n_sites_sel," sites)")))+
      #   theme_bw()+
      #   labs(fill= "Common Vs Rare", title=paste0(species," sampling strategy"), 
      #        y="Alleles (proportion)", x="Number of Sites/ \n Total Samples", colour="")+
      #   geom_hline(yintercept = 0.9, linetype="dashed", alpha=0.5, colour="red")+
      #   geom_hline(yintercept = 0.5, linetype="dashed", alpha=0.5, colour="blue")+
      #   #facet_nested(~t_num_indv+n_sites_sel, scales = "free_x",  space = "free", switch="x")+theme_classic()+
      #   theme(panel.spacing=unit(0,"lines"), axis.title.x = element_blank(),
      #         axis.text = element_text(size=8, face="bold"),
      #         strip.background=element_rect(color="grey30", fill="grey90"), 
      #         panel.grid = element_blank(), legend.position = "bottom",
      #         legend.title = element_blank(),
      #         legend.background = element_rect(fill=alpha("grey", 0.2)))
      # 
      # ggsave(paste0("1. Site Randomisation_Minimum_line.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
      
      
      
      
      
      #OK ready to optimise based off the randomisation!
      
      #optimsie using equal sample numbers based on what was the desired number per site from above
      #then when calcualting the range of alleles captured for each site combination afetr the optimsieation, use all the idnividuals found at a site.
      #this way, you are optimsiing using equal sample numebrs, but then calcualting total allele capture using all samples
      
      
      #idenify how many samples to optimsie for in downstream analyses:
      
      
      auto_nts <- allvals2minsite[which(allvals2minsite$Group=="5% MAF" & allvals2minsite$Group2=="Common" & allvals2minsite$Allele>0.9),] # find the sample combo where the min random allele prop for 5% MAF reaches over 90% common alleles
      auto_nts <- data.frame(auto_nts[order(auto_nts$n_sites_sel),])
      auto_nt_sites <- as.numeric(as.character(auto_nts$n_sites_sel))[1]
      auto_nt_totalsamps <- as.numeric(as.character(auto_nts$t_num_indv ))[1]
      N_t_vec <- auto_nt_sites
      totalsamps <- auto_nt_totalsamps
      sampspersite <- totalsamps/auto_nt_sites
      
    } else {
      N_t_vec <- manual_sites 
      sampspersite <- manual_sampspersite
    }
    
    if (!is.null(sites_to_force) | !is.null(sites_to_exclude)){
      if (any(!length(sites_to_force)<N_t_vec)){
        cat("\n Awwwww, SNAP! you have too many sites to force for the auto_nt option! try using a manual number of sites")
      } else if (any(!length(sites_to_exclude)<N_t_vec)){
        cat("\n Awwwww, SNAP! you have too many sites to exclude for the auto_nt option! try using a manual number of sites")
      }  
    }
    
    poppys <- table(dms$meta$analyses[,site_col_name])
    sampstokeepz <- dms$sample_names[which(dms$meta$analyses[,site_col_name] %in% names(which(poppys>=sampspersite)))]
    pops <- dms$meta$analyses[,site_col_name][which(dms$meta$analyses[,site_col_name] %in% names(which(poppys>=sampspersite)))]
    
    popslatlong <- c()
    for (e in 1:length(unique(pops))){
      popslatlong$lat[e] <-  mean(dms$meta$lat[which(dms$meta$analyses[,site_col_name]==unique(pops)[e])])
      popslatlong$long[e] <-  mean(dms$meta$long[which(dms$meta$analyses[,site_col_name]==unique(pops)[e])])
    }
    popslatlong <- data.frame(popslatlong)
    
    dmssites <- remove.by.list(dms, sampstokeepz)
    dms <- dmssites
    
    
    for (o in 1:length(measurevals)){
      measure <- measurevals[o]
      ulimM <- unlimited_mvals[o]
      allelescapturedfin <- c()
      set.seed(9825)
      sw_out_list <- list()
      
      for ( i in 1:length(N_t_vec) ) {
        N_t <- N_t_vec[i]
        sampspersitesingle <- sampspersite[i]
        
        if (mvalues[o] == "auto"){
          m <- N_t*sampspersitesingle #this sets m to the same value as the optimsied sites*individuals
        } else {
          m <- mvalues[o]
        }
        
        cat("\n Running ", measure," for ", N_t, "Sites  and ",sampspersitesingle ," samps per site...\n")
        
        #now sub sample any sites with more than n samples to nmake optimsiation equal
        samples_df <- dmssites$meta$analyses %>% as.data.frame() %>%
          group_by(!!sym(site_col_name)) %>% slice_sample(n = sampspersitesingle)
        dmssites_temp <- remove.by.list(dmssites,samples_df$sample) #dmssites_temp is now thew dms with equal samples
        gt_sw_comp_sites <- dmssites_temp$gt
        updatedpopsforallelecount <- dmssites_temp$meta$analyses[,site_col_name]
        
        max_wts <- rep(1, length(unique(updatedpopsforallelecount)))
        
        #force or exclude any sites?
        if (!is.null(sites_to_force)&&!is.null(sites_to_exclude)){
          forcedsamps <- which(unique(updatedpopsforallelecount)%in%sites_to_force)
          excludesamps <- which(unique(updatedpopsforallelecount)%in%sites_to_exclude)
          maxws <- replace(max_wts,forcedsamps,0)
          maxws <- replace(maxws,excludesamps,0)
          max_wts[excludesamps] <- 0
          initial_weights <- propose_initial_weights(length(unique(updatedpopsforallelecount)), (N_t-length(sites_to_force)), w_max=maxws)
          initial_weights[forcedsamps] <- 1
          weights_min <- rep(0, length(unique(updatedpopsforallelecount)))
          weights_min <-replace(weights_min,forcedsamps,1)
        } else if (!is.null(sites_to_force)){
          forcedsamps <- which(unique(updatedpopsforallelecount)%in%sites_to_force)
          maxws <- replace(max_wts,forcedsamps,0)
          initial_weights <- propose_initial_weights(length(unique(updatedpopsforallelecount)), (N_t-length(sites_to_force)), w_max=maxws)
          initial_weights[forcedsamps] <- 1
          weights_min <- rep(0, length(unique(updatedpopsforallelecount)))
          weights_min <-replace(weights_min,forcedsamps,1)
        } else if (!is.null(sites_to_exclude)){
          excludesamps <- which(unique(updatedpopsforallelecount)%in%sites_to_exclude)
          max_wts[excludesamps] <- 0
          initial_weights <- propose_initial_weights(length(unique(updatedpopsforallelecount)), N_t, w_max=max_wts)
        } 
        
        
        
        
        if (measure=="psfs"){ 
          #now run the actual psfs optimisation for common alleles
          CommOnly <- gt_sw_comp_sites[,which(colnames(gt_sw_comp_sites)%in%rownames(data.frame(i_sw_common)))]         
          pmac <- gt_to_pop_minor_allele_counts(CommOnly, updatedpopsforallelecount)
          
          opt_results   <- optimize_single_objective( gt=pmac$MAC, N_t=N_t, initial_weights=initial_weights, weights_max=max_wts, 
                                                      measure=measure, max_steps=max_steps, max_t=max_t,
                                                      m=m, p_depends_delta=TRUE, weights_min=weights_min,
                                                      pMAC_mode=TRUE, Nmat=pmac$N, ncpu=1, unlim_m = ulimM, kinall=kinall)
          
        } else {print("uh oh! add more options here from other versions of OptGenMix scripts")}
        
        sw_out_list[[ i ]] <- list(N_t=N_t, m=m, d_opt=opt_results)
        
        OGM_dir_temp <- paste0(OGM_dir,"Temperature_Plots/")
        if (!dir.exists(OGM_dir_temp)) {
          dir.create(OGM_dir_temp, recursive = TRUE)
          cat(paste("Created directory:", OGM_dir_temp, "\n"))
        } else {}
        tiff(paste0(OGM_dir_temp ,N_t,"sites ", sampspersitesingle, "sampspersite ", species, " Temperature Plot T=", max_t, IncludeNA,measure,"m=", m, ".tiff"),
             units = "in", width = 16, height = 10, res = 100)
        par(mfrow = c(1, 1))
        plot(sw_out_list[[i]]$d_opt$value, main= paste0(N_t, "sites ", sampspersitesingle, "sampspersite T_max = ", max_t," ", IncludeNA, " ", measure)) 
        dev.off()
      }
      
      
      
      ##produce table
      solution_table <- mat.or.vec(nrow(pmac$MAC), length(N_t_vec)+3)
      for (i in 1:length(N_t_vec)) { 
        #for each N_t, get the following:
        # 1. get solution table
        solution_table[,i] <- sw_out_list[[i]]$d_opt$weight[max_steps,] 
        
      }
      
      solution_table[,i+1] <- unique(pops)
      solution_table[,i+2] <- popslatlong$lat
      solution_table[,i+3] <- popslatlong$long
      colnames(solution_table) <- c(paste0(N_t_vec,"Si_", sampspersite,"SaSi"), "site", "lat", "long")  #remember to give column names

    }
    
    
    
    
    
    #based on the optimised solution table generated above, calcualte the alelles captured from randomly sampling x samples per site.
    #use the original dms for this, rather than the dms which subsampled sites down to n samples per site
    
    solution_table <- data.frame(solution_table)
    allvals_common2 <- mat.or.vec(max_steps, length(N_t_vec))
    allvals_rare2 <- mat.or.vec(max_steps, length(N_t_vec))
    
    
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
        #SummaryTab$optsitesnsamps[g] <- length(which(solution_table$site==unique(nt_sites)[g] & solution_table[,z]>0)) #identify how many samples for each site 
        SummaryTab$optsitesnsamps[g] <- sampspersite[z]
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
    colnames(allvalsver2) <- c(paste0(N_t_vec, " ", sampspersite), "MAF")
    
    allvalsver22 <- data.frame(melt(allvalsver2, "MAF"))
    colnames(allvalsver22) <- c("MAF","nt","Allele")
    
    allvalsver222 <- separate(allvalsver22, nt, into = c("n_sites_sel", "indv_p_site"), sep = " ")
    
    allvalsver222$t_num_indv <- as.numeric(allvalsver222$n_sites_sel)*as.numeric(allvalsver222$indv_p_site)
    
    write.csv(data.frame(allvalsver222), paste0(OGM_dir, "2. BIGDATA Optimised Site Proportions.csv"),quote=FALSE)
    
    
    #save the range of allele proportion captured when removing samples from optimised combinations
    rm_sample_min2 <- data.frame(allvalsver222 %>% group_by(MAF , n_sites_sel, indv_p_site) %>% slice(which.min(Allele)))
    rm_sample_min2$MinMax <- "Min Value"
    rm_sample_max2 <- data.frame(allvalsver222 %>% group_by(MAF , n_sites_sel, indv_p_site) %>% slice(which.max(Allele)))
    rm_sample_max2$MinMax <- "Max Value"
    rm_sample_range2 <- rbind(rm_sample_min2, rm_sample_max2)
    rm_sample_range2 <- rm_sample_range2[order(rm_sample_range2$MinMax),]
    rm_sample_range2 <- rm_sample_range2[order(rm_sample_range2$MAF),]
    rm_sample_range2 <- rm_sample_range2[order(rm_sample_range2$nt),]
    
    Optvals_extend <- data.frame(matrix(ncol=6, nrow=(nrow(solution_table)-nrow(rm_sample_range2))), stringsAsFactors=FALSE)
    colnames(Optvals_extend) <- colnames(rm_sample_range2)
    Optvals_extend <- rbind(as.data.frame(rm_sample_range2), Optvals_extend)
    Optvals_extend[is.na(Optvals_extend)] <- ""

    solution_tablefin <- cbind(solution_table, Optvals_extend)
    
    write.table(solution_tablefin, paste0(OGM_dir, "2. Optimsised_Site_Combinations.csv"), sep=",",quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    if (auto_nt){
      
      #now plot the variation compared to optimised sample combo - with 5% MAF
      
      GeneralSampComboCommonOnly <- allvalsver22[which(allvalsver22$MAF=="Common"),]
      
      
      findata_common_and_rare$nt <- paste0(N_t_vec, " ", sampspersite)
      comparewithrandomCommon <- findata_common_and_rare[which(findata_common_and_rare$indv_p_site==sampspersite & 
                                                                 findata_common_and_rare$n_sites_sel==N_t_vec &
                                                                 findata_common_and_rare$Group=="5% MAF" &
                                                                 findata_common_and_rare$Group2=="Common"),]
     
      
      minval <- min(c(min(GeneralSampComboCommonOnly$Allele), min(comparewithrandomCommon$Allele))) # find he minimum value to constrain plot
      maxval <- max(c(max(GeneralSampComboCommonOnly$Allele), max(comparewithrandomCommon$Allele))) # find he minimum value to constrain plot
      
      ggplot() +  
        geom_boxplot(data=comparewithrandomCommon, mapping = aes(x = factor(nt), y = Allele, colour="random"), 
                     fill=NA,size=2, position=position_dodge(width=0))+
        #scale_colour_manual(values="black")+
        labs(x = "N Sites", y = "Allele Proportion", colour="Random Sites")+
        #new_scale_color()+
        geom_boxplot(data=GeneralSampComboCommonOnly, mapping = aes(x = factor(nt), y = Allele, group=interaction(nt, MAF), colour="Optimised"), 
                     fill=NA,size=2, position=position_dodge(width=0))+
        scale_colour_manual(values=c("red", "black"))+
        labs(x = "N Sites", y = "Allele Proportion", colour="Optimised Sites")+
        ylim(minval,maxval)+
        theme_minimal()+ 
        ggtitle(paste0("Common AlleleProp Capture when sampling across Optimised sites. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
      
      ggsave(paste0("2. Optmised Site Vs Random Common.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
      
      
      
      #now plot the rare variation compared to optimised sample combo - with 5% MAF
      
      GeneralSampComboRareOnly <- allvalsver22[which(allvalsver22$MAF=="Rare"),]
      
      comparewithrandomRare <- findata_common_and_rare[which(findata_common_and_rare$indv_p_site==sampspersite & 
                                                               findata_common_and_rare$n_sites_sel==N_t_vec &
                                                               findata_common_and_rare$Group=="5% MAF" &
                                                               findata_common_and_rare$Group2=="rare"),]
      
      
      minval <- min(c(min(GeneralSampComboRareOnly$Allele), min(comparewithrandomRare$Allele))) # find he minimum value to constrain plot
      maxval <- max(c(max(GeneralSampComboRareOnly$Allele), max(comparewithrandomRare$Allele))) # find he minimum value to constrain plot
      
      ggplot() +  
        geom_boxplot(data=comparewithrandomRare, mapping = aes(x = factor(nt), y = Allele, colour="Random"), 
                     fill=NA,size=2, position=position_dodge(width=0))+

        labs(x = "N Sites", y = "Allele Proportion", colour="Random Sites")+
        geom_boxplot(data=GeneralSampComboRareOnly, mapping = aes(x = factor(nt), y = Allele, group=interaction(nt, MAF), colour="Optimised"), 
                     fill=NA,size=2, position=position_dodge(width=0))+
        scale_colour_manual(values=c("red", "black"))+
        labs(x = "N Sites", y = "Allele Proportion", colour="Optimised Sites")+
        ylim(minval,maxval)+
        theme_minimal()+ 
        ggtitle(paste0("Common AlleleProp Capture when sampling across Optimised sites. Total Samples: ",nrow(gt_sw_comp), ".Totoal SNPs ",IncludeNA, ": ", (ncol(gt_sw_comp))))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size=20), legend.title = element_text(size=10), legend.text = element_text(size=10), legend.position="right")
      
      ggsave(paste0("2. Optmised Site Vs Random Rare.tiff"), path = paste0(OGM_dir), width = 16, height = 8, dpi = 300, units = "in")
      
    }
  }
}
