# Install requiered packages
# install.packages("stats")
# install.packages("umap")
# install.packages("tidyr")
# install.packages("reshape2")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("hrbrthemes")

# check location library loading information
.libPaths()

args <- commandArgs(TRUE)

# Load required packages
library(stats)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)


path <- args[1]
print("$ ---- Summit based Occupancy matrix  ------- $")
gbin_encode_out <- read.table(paste0(path , args[2]))
gbin_encode_out <- gbin_encode_out %>% distinct()
colnames(gbin_encode_out) <- c("caps", "genomic_site", "presence")
gbin_encode_out_tab <- tidyr::pivot_wider(gbin_encode_out, names_from = caps, values_from = presence)
gbin_encode_out_df <- as.data.frame(gbin_encode_out_tab)

# Replace NA (NA arises when no binding) with zero
gbin_encode_out_df[is.na(gbin_encode_out_df)] = 0
rownames(gbin_encode_out_df) <- gbin_encode_out_df$genomic_site
gbin_encode_out_df <- gbin_encode_out_df[,-1]
gbin_encode_out_mat <- as.matrix(gbin_encode_out_df)
gbin_encode_out_df["boundsites"] <- rowSums(gbin_encode_out_df[1:length(gbin_encode_out_df)])
gbin_encode_out_df <- gbin_encode_out_df[order(-gbin_encode_out_df$boundsites),]

# calculate Percentile of boundsites
quantile(gbin_encode_out_df$boundsites, c(.85, .95)) 

print("Plotting Distribution plot ....") 
pdf(file=paste0(path,"distribution_gbin_encode_out_df.pdf"))

hist(gbin_encode_out_df$boundsites, breaks = 100)

dev.off()

maxbound <- dim(gbin_encode_out_df)[1]
print(paste0("Total sites bound by overall CAPs (in gbinsize bp intervals: ", maxbound))

print("Removing Hotspot regions, using 10,25,50,100 ....")
for (hsrcutoff in c(10,25,50,100)){
  print(hsrcutoff)
  gbin_encode_out_df_hsrfree <- gbin_encode_out_df[which(gbin_encode_out_df$boundsites < hsrcutoff),]
  gbin_encode_out_mat_hsrfree <- as.matrix(gbin_encode_out_df_hsrfree)
  
  print("Counting total sites occupied by each CAP ....")
  total_sites_bound_gbinsummit_hsrfree <- data.frame(total_sites= colSums(gbin_encode_out_df_hsrfree[1:(length(gbin_encode_out_df_hsrfree)-1)]))
  write.table(total_sites_bound_gbinsummit_hsrfree, paste0(path, "total_sites_bound_gbinsummit_hsrfree_",hsrcutoff,".txt"), sep="\t", append = F, quote = F, row.names = F)
  
  
  print("Converting occupancy matrix to pairwise matrix ....")
  converted_matrix_gbinsummit_hsrfree <- crossprod(gbin_encode_out_mat_hsrfree)
  
  write.table(converted_matrix_gbinsummit_hsrfree, paste0(path, "converted_matrix_gbinsummit_hsrfree_",hsrcutoff,".txt"), sep="\t", append = F, quote = F, row.names = F)
  
  # Add overlap and total occupied sites information, so matrix need to be converted to 3 columns data
  
  num_proteins = dim(converted_matrix_gbinsummit_hsrfree)[1]
  unstacked_gbinsummit_hsrfree <- setNames(melt(converted_matrix_gbinsummit_hsrfree), c('CAPb', 'CAPa', 'overlap'))
  total_sites_bound_gbinsummit_hsrfree["CAPa"] <- rownames(total_sites_bound_gbinsummit_hsrfree)
  total_sites_bound_gbinsummit_hsrfree_mat <- merge(total_sites_bound_gbinsummit_hsrfree, unstacked_gbinsummit_hsrfree, by = "CAPa")
  total_sites_bound_gbinsummit_hsrfree_mat_st_pois <- merge(total_sites_bound_gbinsummit_hsrfree, total_sites_bound_gbinsummit_hsrfree_mat, by.x = "CAPa", by.y= "CAPb")
  total_sites_bound_gbinsummit_hsrfree_mat_st_pois <- total_sites_bound_gbinsummit_hsrfree_mat_st_pois[,c(1,3,2,4,5)]
  colnames(total_sites_bound_gbinsummit_hsrfree_mat_st_pois) <- c("CAP_a", "CAP_b", "clean_peak_a", "clean_peak_b", "overlap")
  
  print("Performing Poisson distribution analysis")
  # # calculate coassociation score
  total_sites_bound_gbinsummit_hsrfree_mat_st_pois["new_expected_value"] <- round((as.numeric(total_sites_bound_gbinsummit_hsrfree_mat_st_pois$clean_peak_a) * as.numeric(total_sites_bound_gbinsummit_hsrfree_mat_st_pois$clean_peak_b))/maxbound) #lambda(expected_value)=(n*m)/N
  total_sites_bound_gbinsummit_hsrfree_mat_st_pois["pois_pval"] <- dpois(total_sites_bound_gbinsummit_hsrfree_mat_st_pois$overlap, total_sites_bound_gbinsummit_hsrfree_mat_st_pois$new_expected_value)
  
  #Perform Bonferroni Correction
  total_sites_bound_gbinsummit_hsrfree_mat_st_pois["Corrected_pois_pval"] <- p.adjust(total_sites_bound_gbinsummit_hsrfree_mat_st_pois$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
  
  
  # set the minimum value for p-value which are too low, 323 is lowest value which can be set
  total_sites_bound_gbinsummit_hsrfree_mat_st_pois["Corrected_pois_pval"] <- ifelse(total_sites_bound_gbinsummit_hsrfree_mat_st_pois$Corrected_pois_pval == 0, total_sites_bound_gbinsummit_hsrfree_mat_st_pois$Corrected_pois_pval + 1e-323, total_sites_bound_gbinsummit_hsrfree_mat_st_pois$Corrected_pois_pval)
  total_sites_bound_gbinsummit_hsrfree_mat_st_pois["Coassociation_score"] <- ifelse((total_sites_bound_gbinsummit_hsrfree_mat_st_pois$overlap > total_sites_bound_gbinsummit_hsrfree_mat_st_pois$new_expected_value), -log(total_sites_bound_gbinsummit_hsrfree_mat_st_pois$Corrected_pois_pval,10), log(total_sites_bound_gbinsummit_hsrfree_mat_st_pois$Corrected_pois_pval,10))
  total_sites_bound_gbinsummit_hsrfree_mat_st_pois <- total_sites_bound_gbinsummit_hsrfree_mat_st_pois %>% distinct()
  
  pdf(paste0(path, "histogram_total_sites_bound_gbinsummit_hsrfree_",hsrcutoff,"_mat_st_pois.pdf"))
  hist(total_sites_bound_gbinsummit_hsrfree_mat_st_pois$Coassociation_score, breaks = 100)
  dev.off()
  
  
  # remove below 0 and get only 3 coulmn
  total_sites_bound_gbinsummit_hsrfree_mat_st_pois_re <- total_sites_bound_gbinsummit_hsrfree_mat_st_pois[which(total_sites_bound_gbinsummit_hsrfree_mat_st_pois$Coassociation_score >=0),][,c(1,2,9)]
  colnames(total_sites_bound_gbinsummit_hsrfree_mat_st_pois_re) <- c("CAP_a", "CAP_b", "Score")
  
  print(" Writing output ....")
  write.table(total_sites_bound_gbinsummit_hsrfree_mat_st_pois_re, 
              paste0(args[3],"_",hsrcutoff,"_pois.txt"), 
              sep="\t",
              append = F,
              quote = F,
              row.names = F,
              col.names = F)

}
