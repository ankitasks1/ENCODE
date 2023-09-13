# Load required packages
library(stats)
library(umap)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(hrbrthemes)

# Peak size distribution from all CAPs
MERGED_final_cat_peaking <- fread("/mnt/home3/reid/av638/ENCODE/MERGED_final_cat_peaking.bed")
MERGED_final_cat_peaking <- data.frame(MERGED_final_cat_peaking)
MERGED_final_cat_peaking["Size"] <- MERGED_final_cat_peaking$V3 - MERGED_final_cat_peaking$V2
summary(MERGED_final_cat_peaking$Size)
# Create a histogram
hist(MERGED_final_cat_peaking$Size, freq = FALSE, main = "Histogram and density")

# Calculate density
MERGED_final_cat_peakingSizedx <- density(MERGED_final_cat_peaking$Size)

summary(MERGED_final_cat_peaking$Size)

# Plot the density without histogram
plot(MERGED_final_cat_peakingSizedx, lwd = 2, col = "red",
     main = "Density", xlim = c(0,1000))



ggplot(MERGED_final_cat_peaking, aes(x = Size)) + xlim(c(0,1000)) +
  geom_density(color = "red", # Curve color
               fill = "red",  # Area color
               alpha = 0.5) # Area transparency


# Summit based Occupancy matrix 

#########  400 bp summit data (not peak data). ############
gbin_400_summit_encode_out <- read.table("/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/allouts_1/gbin_400_summit_encode_out.txt")
gbin_400_summit_encode_out <- gbin_400_summit_encode_out %>% distinct()
head(gbin_400_summit_encode_out)
dim(gbin_400_summit_encode_out)
colnames(gbin_400_summit_encode_out) <- c("caps", "genomic_site", "presence")
gbin_400_summit_encode_out_tab <- tidyr::pivot_wider(gbin_400_summit_encode_out, names_from = caps, values_from = presence)
gbin_400_summit_encode_out_df <- as.data.frame(gbin_400_summit_encode_out_tab)
dim(gbin_400_summit_encode_out_df)
# Replace NA (NA arises when no binding) with zero
gbin_400_summit_encode_out_df[is.na(gbin_400_summit_encode_out_df)] = 0
rownames(gbin_400_summit_encode_out_df) <- gbin_400_summit_encode_out_df$genomic_site
gbin_400_summit_encode_out_df <- gbin_400_summit_encode_out_df[,-1]
head(gbin_400_summit_encode_out_df)
dim(gbin_400_summit_encode_out_df)
gbin_400_summit_encode_out_mat <- as.matrix(gbin_400_summit_encode_out_df)
head(gbin_400_summit_encode_out_mat)
gbin_400_summit_encode_out_df["boundsites"] <- rowSums(gbin_400_summit_encode_out_df[1:length(gbin_400_summit_encode_out_df)])
gbin_400_summit_encode_out_df <- gbin_400_summit_encode_out_df[order(-gbin_400_summit_encode_out_df$boundsites),]
dim(gbin_400_summit_encode_out_df)

# calculate Percentile of boundsites
quantile(gbin_400_summit_encode_out_df$boundsites, c(.85, .95)) 
# Distribution plot
ggplot(gbin_400_summit_encode_out_df, aes(x=boundsites)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") + theme_bw() +
  geom_vline(aes(x = boundsites), xintercept = quantile(gbin_400_summit_encode_out_df$boundsites, c(.85, .95)) ,colour=c("green","red"), linetype = "longdash")

ggsave("distribution_gbin_400_summit_encode_out_df.svg", width=12, height=12, units="cm", dpi=96)


#Count total sites occupied by each CAP
total_sites_bound_400_gbin <- data.frame(total_sites= colSums(gbin_400_summit_encode_out_df[1:(length(gbin_400_summit_encode_out_df)-1)]))


# identify DOT: HOT, MILD, COLD
gbin_400_summit_encode_out_df_HOT <- gbin_400_summit_encode_out_df[which(gbin_400_summit_encode_out_df$boundsites >= 30),]
gbin_400_summit_encode_out_df_MILD <- gbin_400_summit_encode_out_df[which(gbin_400_summit_encode_out_df$boundsites >= 8 & gbin_400_summit_encode_out_df$boundsites < 30),]
gbin_400_summit_encode_out_df_COLD <- gbin_400_summit_encode_out_df[which(gbin_400_summit_encode_out_df$boundsites < 8),]

# total sites bound by overall CAPs (in 400 bp intervals)
maxbound <- dim(gbin_400_summit_encode_out_df)[1]

# Remove Hotspot regions, used 0.95 percentile
gbin_400_summit_encode_out_df_hsrfree <- gbin_400_summit_encode_out_df[which(gbin_400_summit_encode_out_df$boundsites < quantile(gbin_400_summit_encode_out_df$boundsites, c(.85, .95))[2]),]
gbin_400_summit_encode_out_mat_hsrfree <- as.matrix(gbin_400_summit_encode_out_df_hsrfree)

#Count total sites occupied by each CAP
total_sites_bound_gbinsummit_400_hsrfree <- data.frame(total_sites= colSums(gbin_400_summit_encode_out_df_hsrfree[1:(length(gbin_400_summit_encode_out_df_hsrfree)-1)]))

write.table(total_sites_bound_gbinsummit_400_hsrfree, "/mnt/beegfs6/home3/reid/av638/ENCODE/total_sites_bound_gbinsummit_400_hsrfree.txt", sep="\t", append = F, quote = F, row.names = F)


# Convert occupancy matrix to pairwise matrix
converted_matrix_gbinsummit_400_hsrfree <- crossprod(gbin_400_summit_encode_out_mat_hsrfree)

write.table(converted_matrix_gbinsummit_400_hsrfree, "/mnt/beegfs6/home3/reid/av638/ENCODE/converted_matrix_gbinsummit_400_hsrfree.txt", sep="\t", append = F, quote = F, row.names = F)


#-------- sum checks for converting -------#
checkloop <- data.frame(gbin_400_summit_encode_out_mat_hsrfree[,c(1,5)])
checkloop["total"] <- rowSums(checkloop)
dim(checkloop[which(checkloop$total == 2),])

# Add overlap and total occupied sites information, so matrix need to be converted to 3 columns data
num_proteins = dim(converted_matrix_gbinsummit_400_hsrfree)[1]
unstacked_gbinsummit_400_hsrfree <- setNames(melt(converted_matrix_gbinsummit_400_hsrfree), c('CAPb', 'CAPa', 'overlap'))
total_sites_bound_gbinsummit_400_hsrfree["CAPa"] <- rownames(total_sites_bound_gbinsummit_400_hsrfree)
total_sites_bound_gbinsummit_400_hsrfree_mat <- merge(total_sites_bound_gbinsummit_400_hsrfree, unstacked_gbinsummit_400_hsrfree, by = "CAPa")
total_sites_bound_gbinsummit_400_hsrfree_mat_st <- merge(total_sites_bound_gbinsummit_400_hsrfree, total_sites_bound_gbinsummit_400_hsrfree_mat, by.x = "CAPa", by.y= "CAPb")
total_sites_bound_gbinsummit_400_hsrfree_mat_st <- total_sites_bound_gbinsummit_400_hsrfree_mat_st[,c(1,3,2,4,5)]
colnames(total_sites_bound_gbinsummit_400_hsrfree_mat_st) <- c("CAP_a", "CAP_b", "clean_peak_a", "clean_peak_b", "overlap")


# sort
total_sites_bound_gbinsummit_400_hsrfree_mat_st <- total_sites_bound_gbinsummit_400_hsrfree_mat_st[order(total_sites_bound_gbinsummit_400_hsrfree_mat_st$CAP_b),]


# # poisson distribution analysis
# # calculate coassociation score
# total_sites_bound_gbinsummit_400_hsrfree_mat_st["new_expected_value"] <- round((as.numeric(total_sites_bound_gbinsummit_400_hsrfree_mat_st$clean_peak_a) * as.numeric(total_sites_bound_gbinsummit_400_hsrfree_mat_st$clean_peak_b))/maxbound) #lambda(expected_value)=(n*m)/N
# total_sites_bound_gbinsummit_400_hsrfree_mat_st["pois_pval"] <- dpois(total_sites_bound_gbinsummit_400_hsrfree_mat_st$overlap, total_sites_bound_gbinsummit_400_hsrfree_mat_st$new_expected_value)
# 
# #Perform Bonferroni Correction
# total_sites_bound_gbinsummit_400_hsrfree_mat_st["Corrected_pois_pval"] <- p.adjust(total_sites_bound_gbinsummit_400_hsrfree_mat_st$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
# total_sites_bound_gbinsummit_400_hsrfree_mat_st["Corrected_pois_pval"] <- ifelse(total_sites_bound_gbinsummit_400_hsrfree_mat_st$Corrected_pois_pval == 0, total_sites_bound_gbinsummit_400_hsrfree_mat_st$Corrected_pois_pval + 1e-323, total_sites_bound_gbinsummit_400_hsrfree_mat_st$Corrected_pois_pval)
# total_sites_bound_gbinsummit_400_hsrfree_mat_st["Coassociation_score"] <- ifelse((total_sites_bound_gbinsummit_400_hsrfree_mat_st$overlap > total_sites_bound_gbinsummit_400_hsrfree_mat_st$new_expected_value), -log(total_sites_bound_gbinsummit_400_hsrfree_mat_st$Corrected_pois_pval,10), log(total_sites_bound_gbinsummit_400_hsrfree_mat_st$Corrected_pois_pval,10))
# total_sites_bound_gbinsummit_400_hsrfree_mat_st <- total_sites_bound_gbinsummit_400_hsrfree_mat_st %>% distinct()
# 
# 
# hist(total_sites_bound_gbinsummit_400_hsrfree_mat_st$Coassociation_score, breaks = 100)
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re <- total_sites_bound_gbinsummit_400_hsrfree_mat_st[,c(1,2,9)]
# colnames(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re) <- c("CAP_a", "CAP_b", "Score")
# 
# # remove POL2 experiments as we removed them from most of the analyses for their different behavior with respect the other TFs
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re <- total_sites_bound_gbinsummit_400_hsrfree_mat_st_re[!str_detect(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re$CAP_a,"^POL"),]
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re <- total_sites_bound_gbinsummit_400_hsrfree_mat_st_re[!str_detect(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re$CAP_b,"^POL"),]
# 
# write.table(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re, 
#             "/mnt/home3/reid/av638/ENCODE/summit_occupancy/total_sites_bound_gbinsummit_400_hsrfree_mat_st_re.txt", 
#             sep="\t",
#             append = F,
#             quote = F,
#             row.names = F,
#             col.names = F)
# 
# extracount <- total_sites_bound_gbinsummit_400_hsrfree_mat_h[which(total_sites_bound_gbinsummit_400_hsrfree_mat_h$Coassociation_score >= 323),]
# extracount1 <- count(extracount, "CAP_a")
# extracount1 <- extracount1[order(-extracount1$freq),]
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_tab <- tidyr::pivot_wider(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re, names_from = CAP_b, values_from = Score) # convert to table
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_df <- as.data.frame(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_tab)
# rownames(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_df) <- total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_df$CAP_a
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat <- as.matrix(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_df[,-1])
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat[total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat>323] <- 323
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat[total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat< (-323)] <- -323
# 
# dim(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat) #498 498
# # Create heatmap
# my_palette <-  colorRampPalette(c("black","whitesmoke", "blue")) (n=650)
# 
# #if variance is zero for some variables (filter-out)
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat <- total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat[apply(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat, 2, sd, na.rm=TRUE) != 0 ,apply(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat, 2, sd, na.rm=TRUE) != 0] 
# 
# dim(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat) # 483 483
# total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_hc <- hclust(as.dist(1-cor(t(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat))), method = "centroid")
# library(gplots)
# heatmap.2(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat, symm = T, Rowv=as.dendrogram(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_hc),
#           dendrogram = "none", col = my_palette,
#           density.info = "none", trace = "none", margins = c(13,0.1), keysize = 0.1,
#           cexRow = 0.9 ,cexCol = 0.9, key.xlab = NA,symbreaks = T,symkey=T,
#           #( "bottom.margin", "left.margin", "top.margin", "right.margin" )
#           key.par=list(mar=c(4.5,0,3,0),cex.main=2,cex.lab=2, cex.axis=2.2),
#           # lmat -- added 2 lattice sections (5 and 6) for padding
#           lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(0.4, 4.8), lwid=c(1, 10, 1))
# 
# 
# # total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat.pdf
# 
# EDASeq::plotPCA(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat)
# 
# breaksListd <- seq(0, 1, by = 0.01)
# library(pheatmap)
# pheatmap_total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_df <- pheatmap(as.matrix(total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_mat), 
#                                                                    na_col = "grey",
#                                                                    color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
#                                                                    clustering_distance_cols = "euclidean",
#                                                                    cluster_rows = T,
#                                                                    cluster_cols = T,
#                                                                    clustering_method = "centroid")
# #save as pheatmap_total_sites_bound_gbinsummit_400_hsrfree_mat_st_re_df.pdf
# 
# 
# # use graph based clustering to identify complexes
# # Cluster data with MCL
# mcxload -abc total_sites_bound_gbinsummit_400_hsrfree_mat_st_re.txt --stream-mirror -write-tab data.tab -o data.mci 
# 
# mcl data.mci -I 1.4
# mcl data.mci -I 2
# mcl data.mci -I 4
# mcl data.mci -I 8
# 
# 
# mcxdump -icl out.data.mci.I14 -tabr data.tab -o dump.data.mci.I14
# mcxdump -icl out.data.mci.I20 -tabr data.tab -o dump.data.mci.I20
# mcxdump -icl out.data.mci.I40 -tabr data.tab -o dump.data.mci.I40
# 
# # Look at cluster sizes
# cat dump.data.mci.I14  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
# cat dump.data.mci.I20  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
# cat dump.data.mci.I40  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
# 
# # Obtain known protein complex data from  http://mips.helmholtz-muenchen.de/corum/#download
# 
# # Filter known complexes to only contain proteins from the ENCODE dataset
# python complex_count_filter.py humanComplexes.txt total_sites_bound_gbinsummit_400_hsrfree_mat_st_re.txt > humanComplexes_filtered.txt
# 
# # Look at accuracy of predicted complexes
# python mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# 
# # Which complexes do we predict perfectly?
# python mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'
# 
# # this script plots the senstivity and specificity to show improvements as inflation parameter increases
# complex_analysis.R
# 
# # http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt
# # cite:  https://www.embopress.org/doi/full/10.15252/msb.20167490
# 
# humap2_complexes_20200809 <- fread("/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809.txt")
# humap2_complexes_20200809 <- data.frame(humap2_complexes_20200809)
# humap2_complexes_20200809["genes"] <- unlist(lapply(strsplit(humap2_complexes_20200809$genenames, " "), function(x) paste0(x, collapse = ",")))
# humap2_complexes_20200809_re <- humap2_complexes_20200809[,c(1,5)]
# write.table(humap2_complexes_20200809_re, "/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809_re.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)
# 
# python benchmarkagain_complex.py dump.data.mci.I40 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# python benchmarkagain_complex.py dump.data.mci.I40 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'
# #1. 1=complex id/name, 2= CAPs present in complex, CAPs identified by MCL, 3=CAPs from 3 intersected with 2, 5= number of known complex members found / total complex members, 6= number of correct complex members / size of predicted complex
# 

# hypergeometric test based analysis (overlap -1 is required as explained here https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/)
total_sites_bound_gbinsummit_400_hsrfree_mat_h <- total_sites_bound_gbinsummit_400_hsrfree_mat_st[,c(1:5)]
total_sites_bound_gbinsummit_400_hsrfree_mat_h["hyper_pval"] <- phyper(total_sites_bound_gbinsummit_400_hsrfree_mat_h$overlap -1, 
                                                                       total_sites_bound_gbinsummit_400_hsrfree_mat_h$clean_peak_a, 
                                                                       maxbound - total_sites_bound_gbinsummit_400_hsrfree_mat_h$clean_peak_a, 
                                                                       total_sites_bound_gbinsummit_400_hsrfree_mat_h$clean_peak_b,
                                                                       lower.tail=FALSE)

#Perform Bonferroni Correction
total_sites_bound_gbinsummit_400_hsrfree_mat_h["Corrected_hyper_pval"] <- p.adjust(total_sites_bound_gbinsummit_400_hsrfree_mat_h$hyper_pval, method = "bonferroni", n = num_proteins * num_proteins) # n = length(pval) or n * n both are same
# set the minimum value for p-value which are too low, 323 is lowest value which can be set
total_sites_bound_gbinsummit_400_hsrfree_mat_h["Corrected_hyper_pval"] <- ifelse(total_sites_bound_gbinsummit_400_hsrfree_mat_h$Corrected_hyper_pval == 0, total_sites_bound_gbinsummit_400_hsrfree_mat_h$Corrected_hyper_pval + 1e-323, total_sites_bound_gbinsummit_400_hsrfree_mat_h$Corrected_hyper_pval)
total_sites_bound_gbinsummit_400_hsrfree_mat_h["Coassociation_score"] <- -log(total_sites_bound_gbinsummit_400_hsrfree_mat_h$Corrected_hyper_pval,10)
total_sites_bound_gbinsummit_400_hsrfree_mat_h <- total_sites_bound_gbinsummit_400_hsrfree_mat_h %>% distinct()

hist(total_sites_bound_gbinsummit_400_hsrfree_mat_h$Coassociation_score, breaks = 100)

ggplot(total_sites_bound_gbinsummit_400_hsrfree_mat_h, aes(x=Coassociation_score)) + 
  geom_histogram( binwidth=5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Histogram of distribution of Coassociation scores") +
  theme_ipsum() +
  theme(plot.title = element_text(size=15))

ggsave("histogram_total_sites_bound_gbinsummit_400_hsrfree_mat_h.svg", width=12, height=12, units="cm", dpi=96)


# get only 3 coulmn
total_sites_bound_gbinsummit_400_hsrfree_mat_h_re <- total_sites_bound_gbinsummit_400_hsrfree_mat_h[,c(1,2,8)]
colnames(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re) <- c("CAP_a", "CAP_b", "Score")

write.table(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re, 
            "/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/total_sites_bound_gbinsummit_400_hsrfree_mat_h_re.txt", 
            sep="\t",
            append = F,
            quote = F,
            row.names = F,
            col.names = F)


total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_tab <- tidyr::pivot_wider(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re, names_from = CAP_b, values_from = Score) # convert to table
total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_df <- as.data.frame(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_tab)
rownames(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_df) <- total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_df$CAP_a
total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat <- as.matrix(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_df[,-1])
total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat[total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat>323] <- 323
total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat[total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat< (-323)] <- -323

dim(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat) #506 506
# Create heatmap
my_palette <-  colorRampPalette(c("black","whitesmoke", "blue")) (n=650)

#if variance is zero for some variables (filter-out)
total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat <- total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat[apply(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat, 2, sd, na.rm=TRUE) != 0 ,apply(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat, 2, sd, na.rm=TRUE) != 0] 

dim(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat) # 491 491
total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_hc <- hclust(as.dist(1-cor(t(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat))), method = "centroid")
library(gplots)
heatmap.2(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat, symm = T, Rowv=as.dendrogram(total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_hc),
          dendrogram = "none", col = my_palette,
          density.info = "none", trace = "none", margins = c(13,0.1), keysize = 0.1,
          cexRow = 0.9 ,cexCol = 0.9, key.xlab = NA,symbreaks = T,symkey=T,
          #( "bottom.margin", "left.margin", "top.margin", "right.margin" )
          key.par=list(mar=c(4.5,0,3,0),cex.main=2,cex.lab=2, cex.axis=2.2),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(0.4, 4.8), lwid=c(1, 10, 1))


# total_sites_bound_gbinsummit_400_hsrfree_mat_h_re_mat.pdf
# 
# # use graph based clustering to identify complexes
# # Cluster data with MCL
# mcxload -abc total_sites_bound_gbinsummit_400_hsrfree_mat_h_re.txt --stream-mirror -write-tab data.tab -o data.mci 
# mcxload -abc gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re.txt --stream-mirror -write-tab data.tab -o data.mci 
# mcl data.mci -I 1.4
# mcl data.mci -I 2
# mcl data.mci -I 4
# mcl data.mci -I 8
# mcl data.mci -I 16
# mcl data.mci -I 20
# mcl data.mci -I 24
# mcl data.mci -I 30
# mcl data.mci -I 40
# mcl data.mci -I 50
# 
# mcxdump -icl out.data.mci.I14 -tabr data.tab -o dump.data.mci.I14
# mcxdump -icl out.data.mci.I20 -tabr data.tab -o dump.data.mci.I20
# mcxdump -icl out.data.mci.I40 -tabr data.tab -o dump.data.mci.I40
# mcxdump -icl out.data.mci.I80 -tabr data.tab -o dump.data.mci.I80
# mcxdump -icl out.data.mci.I160 -tabr data.tab -o dump.data.mci.I160
# mcxdump -icl out.data.mci.I200 -tabr data.tab -o dump.data.mci.I200
# 
# # Look at cluster sizes
# cat dump.data.mci.I14  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
# cat dump.data.mci.I20  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
# cat dump.data.mci.I40  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
# 
# # details of what each python scripts do are here:
# # Filter known complexes to only contain proteins from the ENCODE dataset
# python3 complex_count_filter.py humanComplexes.txt total_sites_bound_gbinsummit_400_hsrfree_mat_h_re.txt > humanComplexes_filtered.txt
# 
# # Look at accuracy of predicted complexes
# python3 mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# 
# # Which complexes do we predict perfectly?
# python3 mod_benchmark_complexes.py dump.data.mci.I160 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'
# 
# # this script plots the senstivity and specificity to show improvements as inflation parameter increases
# complex_analysis.R
# 
# # http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt
# # cite:  https://www.embopress.org/doi/full/10.15252/msb.20167490
# 
# humap2_complexes_20200809 <- fread("/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809.txt")
# humap2_complexes_20200809 <- data.frame(humap2_complexes_20200809)
# humap2_complexes_20200809["genes"] <- unlist(lapply(strsplit(humap2_complexes_20200809$genenames, " "), function(x) paste0(x, collapse = ",")))
# humap2_complexes_20200809_re <- humap2_complexes_20200809[,c(1,5)]
# write.table(humap2_complexes_20200809_re, "/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809_re.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)
# 
# python3 benchmarkagain_complex.py dump.data.mci.I160 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# python3 benchmarkagain_complex.py dump.data.mci.I160 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'
# #1. 1=complex id/name, 2= CAPs present in complex, CAPs identified by MCL, 3=CAPs from 3 intersected with 2, 5= number of known complex members found / total complex members, 6= number of correct complex members / size of predicted complex
# 
# python3 ./scripts/mod_benchmark_complexes.py dump.data.mci.I160 humanComplexes_filtered.txt > interactions_CAPs_K562_corum.txt
# python3 ./scripts/mod_benchmark_complexes.py dump.data.mci.I160 humap2_complexes_20200809_re.txt  > interactions_CAPs_K562_humap2.txt

# running inside /mnt/home3/reid/av638/ENCODE/summit_occupancy/K562




# Obtain known protein complex data from  http://mips.helmholtz-muenchen.de/corum/#download
# overlap with databases
#corum
python3 ./scripts/corum_complex_count_filter.py humanComplexes.txt gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re.txt > humanComplexes_filtered.txt

#humap2
python3 ./scripts/humap2_complex_count_filter.py humap2_complexes_20200809.txt gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re.txt  > humap2_complexes_20200809_filtered.txt

#both corum and humap2 complex benchmarking in single script 
python3 ./scripts/mcl_benchmark.py gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re.txt


# k562, corum
i_val_k562_corum = list()
for (intfiles in list.files(path = "/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562",pattern = "K562_corum", full.names = T)){
  print(intfiles)
  temp <- fread(intfiles)
  colnames(temp) <- c("id", "database", "mclprediction", "interaction", "sensitivity", "specificity")
  # print(head(temp))
  temp <- data.frame(temp[,c(1,5,6)])
  temp["I_val"] <- as.numeric(i)/10
  i = gsub(".txt","",gsub("/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/interactions.mci.CAPs_K562_corum_","",intfiles))
  print(i)
  i_val_k562_corum[[paste0('I',as.numeric(i)/10)]] <- temp
}
i_val_k562_corum_df <- do.call(rbind.data.frame, i_val_k562_corum)
plot(i_val_k562_corum_df$I_val, i_val_k562_corum_df$sensitivity)

# k562, humap2

i_val_k562_humap2 = list()
for (intfiles in list.files(path = "/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562",pattern = "K562_humap2", full.names = T)){
  print(intfiles)
  temp <- fread(intfiles)
  colnames(temp) <- c("id", "database", "mclprediction", "interaction", "sensitivity", "specificity")
  # print(head(temp))
  temp <- data.frame(temp[,c(1,5,6)])
  temp["I_val"] <- as.numeric(i)/10
  i = gsub(".txt","",gsub("/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/interactions.mci.CAPs_K562_humap2_","",intfiles))
  print(i)
  i_val_k562_humap2[[paste0('I',as.numeric(i)/10)]] <- temp
}
i_val_k562_humap2_df <- do.call(rbind.data.frame, i_val_k562_humap2)
plot(i_val_k562_humap2_df$I_val, i_val_k562_humap2_df$sensitivity)


# pairs benchmarking
#stringDB
# filter cutoff is required to determine the interacting partner, I use > 900 score
library(data.table)
protein9606_id_links <- fread("/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/9606.protein.links.v12.0.txt")
protein9606_id_links <- data.frame(protein9606_id_links[order(-protein9606_id_links$combined_score),])
hist(protein9606_id_links$combined_score, breaks = 100)
dim(protein9606_id_links)#choose cutoff 900, why? highly interacting pairs
dim(protein9606_id_links[which(protein9606_id_links$combined_score > 900),])
protein9606_id_links_high <- protein9606_id_links[which(protein9606_id_links$combined_score > 900),]
write.table(protein9606_id_links_high,
            "/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/9606.protein_id_links_high.txt",
            sep = " ",
            append = F,
            quote = F,
            col.names = T,
            row.names = F)

# take different cutoffs of stringDB scores and writeout
for (exp_scores in seq(200,1000,100)){
  print(exp_scores)
  temp <- protein9606_id_links[which(protein9606_id_links$combined_score > exp_scores),]
  write.table(temp,
              paste0("/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/9606.protein_id_links_",exp_scores,".txt"),
                     sep = " ",
                     append = F,
                     quote = F,
                     col.names = T,
                     row.names = F)

}


# filter interaction of encode data based on cutoff
gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re <- fread("/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re.txt")
gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re <- data.frame(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re)
summary(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re$V3)
hist(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re$V3, breaks = 100)
gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re["id"] <- paste0(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re$V1, "_" ,
                                                                      gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re$V2)

# get only those CAPs which are in encode data
python3 ./scripts/stringDB_complex_count_filter.py 9606.protein.info.v12.0.txt 9606.protein.links.v12.0.txt gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re.txt > protein_id_links_filtered.txt
python3 ./scripts/stringDB_complex_count_filter.py 9606.protein.info.v12.0.txt 9606.protein_id_links_high.txt gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re.txt > protein_id_links_high_filtered.txt
# for multiple cutoff of strdb
./strdb.sh
./strdb_encode.sh 

strdb_encode_list <- list()
for (files in list.files(path="/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562",pattern = "pr.txt")){
  print(files)
  filescore = gsub("_filtered_pr.txt","",gsub("protein_id_links_","",files))
  temp <- read.table(paste0("/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/",files))
  temp["score"] <- filescore
  strdb_encode_list[[filescore]] <- temp
}

strdb_encode_list_df <- data.frame(do.call(rbind.data.frame, strdb_encode_list))

  
ggplot(strdb_encode_list_df,aes(x=V3,y=V2,col=score, label=score))+
    geom_line(aes(color=score))+
    geom_point(aes(color=score), shape=20)+
    scale_color_manual(values=c("darkgrey","black","blue","red","green","orange","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
    theme_bw()
  
  

# assign id to stringDB all proteins interactions
python3 ./scripts/stringDB_complex_id_assign.py 9606.protein.info.v12.0.txt 9606.protein.links.v12.0.txt 

protein_id_links_id_assigned <- read.table("/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/protein_id_links_id_assigned.txt")
protein_id_links_id_assigned["id"] <- paste0(protein_id_links_id_assigned$V1, "_", 
                                             protein_id_links_id_assigned$V2)

# those which are present in encode data
protein_id_links_id_filtered <- read.table("/mnt/home3/reid/av638/ENCODE/summit_occupancy/K562/protein_id_links_filtered.txt")
protein_id_links_id_filtered_assigned <- cSplit(protein_id_links_id_filtered, "V2",",")
protein_id_links_id_filtered_assigned <- data.frame(protein_id_links_id_filtered_assigned)
protein_id_links_id_filtered_assigned["id"] <- paste0(protein_id_links_id_filtered_assigned$V2_1, "_", 
                                                      protein_id_links_id_filtered_assigned$V2_2)
protein_id_links_id_filtered_assigned <- protein_id_links_id_filtered_assigned[,c(3,4,2,5)]

gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB <- merge(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re, protein_id_links_id_assigned, by="id", all.x = TRUE)
# remove the both side same CAPs
gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB <- gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB[which(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB$V1.x != gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB$V2.x),]
colnames(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB) <- c("id",
                                                                           "CAP_a",
                                                                           "CAP_b",
                                                                           "Score",
                                                                           "Protein1",
                                                                           "Protein2",
                                                                           "StringDBscore")


# To understand predicted and benchmark pairs, a threshold is required for both predicted-coassosation score and benchmark-stringDB score
gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_int <-
  gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re[which(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_re$V3 > 200),]

protein_id_links_id_filtered_assigned_int <- protein_id_links_id_filtered_assigned[which(protein_id_links_id_filtered_assigned$V3 > 500),]

gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB_int <- merge(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_int, protein_id_links_id_filtered_assigned_int, by="id", all.x = TRUE)
# remove a _ a
# remove a vs b and b vs a

colnames(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB_int) <- c("id",
                                                                           "CAP_a",
                                                                           "CAP_b",
                                                                           "Score",
                                                                           "Protein1",
                                                                           "Protein2",
                                                                           "StringDBscore")


predictaed_pairs_cs200 <- dim(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_int)[1]
behcmark_pairs_stringDB_s500 <- dim(protein_id_links_id_filtered_assigned_int)[1]
overlap_predicted_benchmark_stringDB <- dim(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB_int[!is.na(gbin_400_total_sites_bound_summit_hsrfree_mat_st_h_stringDB_int$Protein1),])[1]
# 
# TP= overlap between encode and stringDB
# FN= only present in stringDB
# FP=only in encode
# TN=neither in encode nore in stringDB, all other non interactions 
# 
# calculate: PR-curve: higher precision and recall and F-score
#------------  END OF ANALYSIS -----------#


### peak based 1000bp occupancy matrix

gbin_1000_peaks_encode_out <- read.table("/mnt/home3/reid/av638/ENCODE/pairwise/gbin_1000_peaks_encode_out.txt")
gbin_1000_peaks_encode_out <- gbin_1000_peaks_encode_out %>% distinct()
head(gbin_1000_peaks_encode_out)
dim(gbin_1000_peaks_encode_out)
colnames(gbin_1000_peaks_encode_out) <- c("caps", "genomic_site", "presence")
gbin_1000_peaks_encode_out_tab <- tidyr::pivot_wider(gbin_1000_peaks_encode_out, names_from = caps, values_from = presence)
gbin_1000_peaks_encode_out_df <- as.data.frame(gbin_1000_peaks_encode_out_tab)

# Replace NA (NA arises when no binding) with zero
gbin_1000_peaks_encode_out_df[is.na(gbin_1000_peaks_encode_out_df)] = 0
rownames(gbin_1000_peaks_encode_out_df) <- gbin_1000_peaks_encode_out_df$genomic_site
gbin_1000_peaks_encode_out_df <- gbin_1000_peaks_encode_out_df[,-1]
head(gbin_1000_peaks_encode_out_df)
dim(gbin_1000_peaks_encode_out_df)
gbin_1000_peaks_encode_out_mat <- as.matrix(gbin_1000_peaks_encode_out_df)
head(gbin_1000_peaks_encode_out_mat)
gbin_1000_peaks_encode_out_df["boundsites"] <- rowSums(gbin_1000_peaks_encode_out_df[1:length(gbin_1000_peaks_encode_out_df)])
gbin_1000_peaks_encode_out_df <- gbin_1000_peaks_encode_out_df[order(-gbin_1000_peaks_encode_out_df$boundsites),]
dim(gbin_1000_peaks_encode_out_df)


gbin_1000_peaks_encode_out_df_hsrfree <- gbin_1000_peaks_encode_out_df[which(gbin_1000_peaks_encode_out_df$boundsites <= 50),]
gbin_1000_peaks_encode_out_df_hsrfree["coordinates"] <- rownames(gbin_1000_peaks_encode_out_df_hsrfree)

write.table(data.frame(gbin_1000_peaks_encode_out_df_hsrfree$coordinates),
            "gbin_1000_peaks_encode_out_df_hsrfree.txt",
            sep ="\t",
            quote=F,
            append=F,
            row.names = F,
            col.names = F)
#Count total sites occupied
total_sites_bound_gbin_1000 <- data.frame(total_sites= colSums(gbin_1000_peaks_encode_out_df[1:(length(gbin_1000_peaks_encode_out_df)-1)]))

# Distribution plot
ggplot(gbin_1000_peaks_encode_out_df, aes(x=boundsites)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") + theme_bw() +
  geom_vline(aes(x = boundsites), xintercept = quantile(gbin_1000_peaks_encode_out_df$boundsites, c(.85, .95)) ,colour=c("green","red"), linetype = "longdash")

# calculate Percentile of boundsites
quantile(gbin_1000_peaks_encode_out_df$boundsites, c(.85, .95)) 

# identify DOT: HOT, MILD, COLD
gbin_1000_peaks_encode_out_df_HOT <- gbin_1000_peaks_encode_out_df[which(gbin_1000_peaks_encode_out_df$boundsites >= 47),]
gbin_1000_peaks_encode_out_df_MILD <- gbin_1000_peaks_encode_out_df[which(gbin_1000_peaks_encode_out_df$boundsites >= 11 & gbin_1000_peaks_encode_out_df$boundsites < 47),]
gbin_1000_peaks_encode_out_df_COLD <- gbin_1000_peaks_encode_out_df[which(gbin_1000_peaks_encode_out_df$boundsites < 11),]

# get count distribution data out
count_gbin_1000_peaks_encode_out_df <- count(gbin_1000_peaks_encode_out_df$boundsites)
plot(count_gbin_1000_peaks_encode_out_df$x, count_gbin_1000_peaks_encode_out_df$freq)

# prepare hierarchical cluster
gbin_1000_peaks_encode_out_hc = hclust(dist(t(gbin_1000_peaks_encode_out_df[,c(1:(length(gbin_1000_peaks_encode_out_df) -1))])))
plot(gbin_1000_peaks_encode_out_hc)
# gbin_1000_peaks_encode_out_hc.pdf

# separate groups
gbin_1000_peaks_encode_out_cutree <- data.frame(groups=cutree(gbin_1000_peaks_encode_out_hc, k = 50))
gbin_1000_peaks_encode_out_cutree["caps"] <- rownames(gbin_1000_peaks_encode_out_cutree)
count_gbin_1000_peaks_encode_out_cutree <- count(gbin_1000_peaks_encode_out_cutree$groups)
colnames(count_gbin_1000_peaks_encode_out_cutree) <- c("complex", "num_of_caps")
plot(count_gbin_1000_peaks_encode_out_cutree$complex, count_gbin_1000_peaks_encode_out_cutree$num_of_caps)

gbin_1000_peaks_encode_out_cutree_merge <- ddply(gbin_1000_peaks_encode_out_cutree, .(groups), summarize, caps = toString(caps))
# gbin_1000_peaks_encode_out_cutree_merge[gbin_1000_peaks_encode_out_cutree_merge$caps %like% "CTCF", ]
gbin_1000_peaks_encode_out_cutree_merge[grepl("\\CTCF\\b", gbin_1000_peaks_encode_out_cutree_merge$caps, ignore.case = TRUE),]

# Convert occupancy matrix to pairwise matrix
# Create empty matrix to store data
converted_matrix_gbin_1000_peaks <- matrix(0, nrow = ncol(gbin_1000_peaks_encode_out_mat), ncol = ncol(gbin_1000_peaks_encode_out_mat))
rownames(converted_matrix_gbin_1000_peaks)  <- colnames(gbin_1000_peaks_encode_out_mat)
colnames(converted_matrix_gbin_1000_peaks)  <- colnames(gbin_1000_peaks_encode_out_mat)

# Run loop to create pairwise matrix
for (i in 1:ncol(gbin_1000_peaks_encode_out_mat)){
  for  (j in 1:ncol(gbin_1000_peaks_encode_out_mat)){
    total_sites <- sum(gbin_1000_peaks_encode_out_mat[,i])
    print(colnames(converted_matrix_gbin_1000_peaks[,c(i,j)]))
    # print(total_sites)
    # print(paste0(i,j))
    # print(paste0(gbin_1000_peaks_encode_out_mat[,i],gbin_1000_peaks_encode_out_mat[,j]))
    # print(paste0(sum(as.logical(ifelse(gbin_1000_peaks_encode_out_mat[,i] == 1 & gbin_1000_peaks_encode_out_mat[,j] == 1, 1,0))), '/',total_sites))
    converted_matrix_gbin_1000_peaks[i,j] <- sum(as.logical(ifelse(gbin_1000_peaks_encode_out_mat[,i] == 1 & gbin_1000_peaks_encode_out_mat[,j] == 1, 1,0)))
  }
}

write.table(converted_matrix_gbin_1000_peaks, "/mnt/beegfs6/home3/reid/av638/ENCODE/pairwise/converted_matrix_gbin_1000_peaks.txt", sep="\t", append = F, quote = F, row.names = F)
# 
breaksListd <- seq(0, 100, by = 0.01)
library(pheatmap)
pheatmap_converted_matrix_gbin_1000_peaks <- pheatmap(as.matrix(converted_matrix_gbin_1000_peaks),
                                                              na_col = "grey",
                                                              color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
                                                              clustering_distance_cols = "euclidean",
                                                              cluster_rows = T,
                                                              cluster_cols = T,
                                                              clustering_method = "ward.D")
#save as pheatmap_converted_matrix_gbin_1000_peaks.pdf

# prepare hierarchical cluster
converted_matrix_gbin_1000_peaks_hc = hclust(dist(t(converted_matrix_gbin_1000_peaks)))
plot(converted_matrix_gbin_1000_peaks_hc)
# converted_matrix_gbin_1000_peaks_hc.pdf
# Add overlap and total occupied sites information
library(reshape2)
library(dplyr)
num_proteins = 505
min_peaks = 100
unstacked_gbin_1000_peaks <- setNames(melt(converted_matrix_gbin_1000_peaks), c('CAPb', 'CAPa', 'overlap'))
total_sites_bound_gbin_1000_re <- total_sites_bound_gbin_1000
total_sites_bound_gbin_1000_re["CAPa"] <- rownames(total_sites_bound_gbin_1000_re)
total_sites_bound_gbin_1000_mat <- merge(total_sites_bound_gbin_1000_re, unstacked_gbin_1000_peaks, by = "CAPa")
total_sites_bound_gbin_1000_mat_st <- merge(total_sites_bound_gbin_1000_re, total_sites_bound_gbin_1000_mat, by.x = "CAPa", by.y= "CAPb")
total_sites_bound_gbin_1000_mat_st <- total_sites_bound_gbin_1000_mat_st[,c(1,3,2,4,5)]
colnames(total_sites_bound_gbin_1000_mat_st) <- c("CAP_a", "CAP_b", "clean_peak_a", "clean_peak_b", "overlap")

# total sites bound by overall CAPs
maxbound <- dim(gbin_1000_peaks_encode_out_mat)[1]

total_sites_bound_gbin_1000_mat_st["new_expected_value"] <- round((as.numeric(total_sites_bound_gbin_1000_mat_st$clean_peak_a) * as.numeric(total_sites_bound_gbin_1000_mat_st$clean_peak_b))/maxbound) #lambda(expected_value)=(n*m)/N
total_sites_bound_gbin_1000_mat_st["pois_pval"] <- dpois(total_sites_bound_gbin_1000_mat_st$overlap, total_sites_bound_gbin_1000_mat_st$new_expected_value)
#Perform Bonferroni Correction
total_sites_bound_gbin_1000_mat_st["Corrected_pois_pval"] <- p.adjust(total_sites_bound_gbin_1000_mat_st$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
total_sites_bound_gbin_1000_mat_st["Corrected_pois_pval"] <- ifelse(total_sites_bound_gbin_1000_mat_st$Corrected_pois_pval == 0, total_sites_bound_gbin_1000_mat_st$Corrected_pois_pval + 1e-300, total_sites_bound_gbin_1000_mat_st$Corrected_pois_pval)
total_sites_bound_gbin_1000_mat_st["Coassociation_score"] <- ifelse((total_sites_bound_gbin_1000_mat_st$overlap > total_sites_bound_gbin_1000_mat_st$new_expected_value), -log(total_sites_bound_gbin_1000_mat_st$Corrected_pois_pval,10), log(total_sites_bound_gbin_1000_mat_st$Corrected_pois_pval,10))
total_sites_bound_gbin_1000_mat_st <- total_sites_bound_gbin_1000_mat_st %>% distinct()


hist(total_sites_bound_gbin_1000_mat_st$Coassociation_score, breaks = 100)

#filter low count peaks
#sideA
total_sites_bound_gbin_1000_mat_st_filt <- total_sites_bound_gbin_1000_mat_st[which(total_sites_bound_gbin_1000_mat_st$clean_peak_a > min_peaks),]
#sideB
total_sites_bound_gbin_1000_mat_st_filt <- total_sites_bound_gbin_1000_mat_st_filt[which(total_sites_bound_gbin_1000_mat_st_filt$clean_peak_b > min_peaks),]
# 
# #remove self-hits same protein on both sides
# total_sites_bound_gbin_1000_mat_st_filt <- total_sites_bound_gbin_1000_mat_st_filt[which(total_sites_bound_gbin_1000_mat_st_filt$CAP_a != total_sites_bound_gbin_1000_mat_st_filt$CAP_b),]
# 
# #remove 0s coassociation score
# total_sites_bound_gbin_1000_mat_st_filt2 <- total_sites_bound_gbin_1000_mat_st_filt[which(total_sites_bound_gbin_1000_mat_st_filt$Coassociation_score != 0),]
# 
# #remove +/- 1.4 coassociation
# total_sites_bound_gbin_1000_mat_st_filt3 <- total_sites_bound_gbin_1000_mat_st_filt2[which(total_sites_bound_gbin_1000_mat_st_filt2$Coassociation_score > 1.4 | total_sites_bound_gbin_1000_mat_st_filt2$Coassociation_score < -1.4),]
# write.table(total_sites_bound_gbin_1000_mat_st_filt3, "/mnt/beegfs6/home3/reid/av638/ENCODE/pairwise/total_sites_bound_gbin_1000_mat_st_filt3.txt", sep="\t", append = F, quote = F, row.names = F)

total_sites_bound_gbin_1000_mat_st_re <- total_sites_bound_gbin_1000_mat_st_filt[,c(1,2,9)]

# total_sites_bound_gbin_1000_mat_st_re <- total_sites_bound_gbin_1000_mat_st_filt3[,c(1,2,9)]
colnames(total_sites_bound_gbin_1000_mat_st_re) <- c("CAP_a", "CAP_b", "score")
write.table(total_sites_bound_gbin_1000_mat_st_re, "/mnt/beegfs6/home3/reid/av638/ENCODE/pairwise/total_sites_bound_gbin_1000_mat_st_re.txt", sep="\t", append = F, quote = F, row.names = F)
hist(total_sites_bound_gbin_1000_mat_st_re$score, breaks = 100)

total_sites_bound_gbin_1000_mat_st_re_tab <- tidyr::pivot_wider(total_sites_bound_gbin_1000_mat_st_re, names_from = CAP_b, values_from = score)

#NA will come for those with no combinations matched
total_sites_bound_gbin_1000_mat_st_re_df <- as.data.frame(total_sites_bound_gbin_1000_mat_st_re_tab)
rownames(total_sites_bound_gbin_1000_mat_st_re_df) <- total_sites_bound_gbin_1000_mat_st_re_df$CAP_a
total_sites_bound_gbin_1000_mat_st_re_df <- total_sites_bound_gbin_1000_mat_st_re_df[,-1]
dim(total_sites_bound_gbin_1000_mat_st_re_df)

#replace NA with 0 to show the data was obsent
total_sites_bound_gbin_1000_mat_st_re_df[is.na(total_sites_bound_gbin_1000_mat_st_re_df)] <- 0

dim(total_sites_bound_gbin_1000_mat_st_re_df)

EDASeq::plotPCA(as.matrix(total_sites_bound_gbin_1000_mat_st_re_df))

breaksListd <- seq(0, 1, by = 0.01)
library(pheatmap)
pheatmap_total_sites_bound_gbin_1000_mat_st_re_df <- pheatmap(as.matrix(total_sites_bound_gbin_1000_mat_st_re_df), 
                                               na_col = "grey",
                                               color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
                                               clustering_distance_cols = "euclidean",
                                               cluster_rows = T,
                                               cluster_cols = T,
                                               clustering_method = "ward.D")
#save as pheatmap_total_sites_bound_gbin_1000_mat_st_re_df.pdf

# prepare hierarchical cluster
total_sites_bound_gbin_1000_mat_st_re_hc = hclust(dist(t(total_sites_bound_gbin_1000_mat_st_re_df)))
plot(total_sites_bound_gbin_1000_mat_st_re_hc)
# total_sites_bound_gbin_1000_mat_st_re_hc.pdf

# gbin_1000_peaks_encode_out_df_suz12 <- gbin_1000_peaks_encode_out_df[,c("SUZ12", "EZH2", "CBX2", "CBX8", "SIX5")]
# EDASeq::plotPCA(as.matrix(gbin_1000_peaks_encode_out_df_suz12))
# # test pattern of bidning
# plot(hclust(dist(t(gbin_1000_peaks_encode_out_df_suz12[,c(1:length(gbin_1000_peaks_encode_out_df_suz12))]))))

# Assuming your data is stored in a matrix called 'data_matrix'
# Perform PCA
#Heatmap_encode_ivalab_countout_df_filt.pdf
EDASeq::plotPCA(as.matrix(gbin_1000_peaks_encode_out_mat))
#transpose
tgbin_1000_peaks_encode_out_mat <- t(gbin_1000_peaks_encode_out_mat)
dim(tgbin_1000_peaks_encode_out_mat)
rownames(tgbin_1000_peaks_encode_out_mat)
tgbin_1000_peaks_encode_out_mat = data.frame(tgbin_1000_peaks_encode_out_mat)
tgbin_1000_peaks_encode_out_mat["Color"] <- rownames(tgbin_1000_peaks_encode_out_mat)
dim(tgbin_1000_peaks_encode_out_mat)
tgbin_1000_peaks_encode_out_matdfx <- tgbin_1000_peaks_encode_out_mat[c(1:(length(tgbin_1000_peaks_encode_out_mat) -1))]
gbin_1000_peaks_encod_pca_result <- prcomp(tgbin_1000_peaks_encode_out_matdfx, center = TRUE, scale. = TRUE)
#PC<-prcomp(dfx, scale. = T)
gbin_1000_peaks_encod_pca_PCi <- data.frame(gbin_1000_peaks_encod_pca_result$x, Color=tgbin_1000_peaks_encode_out_mat$Color)
percentage <- round(gbin_1000_peaks_encod_pca_result$sdev^2 / sum(gbin_1000_peaks_encod_pca_result$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentage <- paste(colnames(gbin_1000_peaks_encod_pca_PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

ggplot(gbin_1000_peaks_encod_pca_PCi,aes(x=PC1,y=PC2, color=Color, label=Color))+
   xlab(percentage[1]) + ylab(percentage[2])+
  geom_point(size=1,alpha=1) +theme_bw() + xlim(c(-100,500)) + ylim(c(-500,100))+
  geom_text(size=1,label=gbin_1000_peaks_encod_pca_PCi$Color, check_overlap = T)


tgbin_1000_peaks_encode_out_matdfx_umap_data <- umap(tgbin_1000_peaks_encode_out_matdfx)
tgbin_1000_peaks_encode_out_matdfx_umap_df <- data.frame(tgbin_1000_peaks_encode_out_matdfx_umap_data$layout)
colnames(tgbin_1000_peaks_encode_out_matdfx_umap_df) <- c("V1", "V2")
tgbin_1000_peaks_encode_out_matdfx_umap_df["Color"] <- rownames(tgbin_1000_peaks_encode_out_matdfx_umap_df)

ggplot(tgbin_1000_peaks_encode_out_matdfx_umap_df,aes(x=V1,y=V2, label=Color))+
  #geom_point(size=0.1,alpha=1) +
  theme_bw() +
  geom_text(size=4,label=tgbin_1000_peaks_encode_out_matdfx_umap_df$Color, check_overlap = T)
#tgbin_1000_peaks_encode_out_matdfx_umap_df.pdf

# Assuming you have 'k' principal components stored in 'pca_result$u'
gbin_1000_peaks_encod_umap_result <- umap(gbin_1000_peaks_encod_pca_result$x[, 1:20])
colnames(gbin_1000_peaks_encod_umap_result$layout) <- c("UMAP1", "UMAP2")
gbin_1000_peaks_encod_umap_result_layout <- data.frame(gbin_1000_peaks_encod_umap_result$layout)
gbin_1000_peaks_encod_umap_result_layout["Color"] <- rownames(gbin_1000_peaks_encod_umap_result_layout)

# Plot UMAP
plot(gbin_1000_peaks_encod_umap_result$layout[, 1], gbin_1000_peaks_encod_umap_result$layout[, 2], type = "n", xlab = "UMAP 1", ylab = "UMAP 2")
points(gbin_1000_peaks_encod_umap_result$layout[, 1], gbin_1000_peaks_encod_umap_result$layout[, 2], col = "red")

ggplot(gbin_1000_peaks_encod_umap_result_layout,aes(x=UMAP1,y=UMAP2, label=Color))+
  #geom_point(size=0.1,alpha=1) +
  theme_bw() +
  geom_text(size=0.4,label=gbin_1000_peaks_encod_umap_result_layout$Color, check_overlap = T)

#umap_gbin_1000_peaks_encod_umap_result_layout.svg
#umap_gbin_1000_peaks_encod_umap_result_layout.pdf

# pairwise, overlap matrix

num_proteins = 794
intersected_one_all <- fread("/mnt/beegfs6/home3/reid/av638/ENCODE/pairwise/merged_out/intersected_one_all.txt")
intersected_one_all <- data.frame(intersected_one_all)
colnames(intersected_one_all) <- c("CAP_a", "CAP_b", "peaks_a", "peaks_b", "clean_peak_a", "clean_peak_b","overlap")
intersected_one_all["new_expected_value"] <- round((as.numeric(intersected_one_all$clean_peak_a) * as.numeric(intersected_one_all$clean_peak_b))/250000) #lambda(expected_value)=(n*m)/N
intersected_one_all["pois_pval"] <- dpois(intersected_one_all$overlap, intersected_one_all$new_expected_value)
#Perform Bonferroni Correction
intersected_one_all["Corrected_pois_pval"] <- p.adjust(intersected_one_all$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
intersected_one_all["Corrected_pois_pval"] <- ifelse(intersected_one_all$Corrected_pois_pval == 0, intersected_one_all$Corrected_pois_pval + 1e-300, intersected_one_all$Corrected_pois_pval)
intersected_one_all["Coassociation_score"] <- ifelse((intersected_one_all$overlap > intersected_one_all$new_expected_value), -log(intersected_one_all$Corrected_pois_pval,10), log(intersected_one_all$Corrected_pois_pval,10))
intersected_one_all <- intersected_one_all %>% distinct()
head(intersected_one_all)
#filter low count peaks
#sideA
intersected_one_all_filt <- intersected_one_all[which(intersected_one_all$clean_peak_a > 100),]
#sideB
intersected_one_all_filt <- intersected_one_all_filt[which(intersected_one_all_filt$clean_peak_b > 100),]

#remove self-hits same protein on both sides
intersected_one_all_filt <- intersected_one_all_filt[which(intersected_one_all_filt$CAP_a != intersected_one_all_filt$CAP_b),]

#remove 0s coassociation score
intersected_one_all_filt2 <- intersected_one_all_filt[which(intersected_one_all_filt$Coassociation_score != 0),]

#remove +/- 1.4 coassociation
intersected_one_all_filt3 <- intersected_one_all_filt2[which(intersected_one_all_filt2$Coassociation_score > 1.4 | intersected_one_all_filt2$Coassociation_score < -1.4),]
write.table(intersected_one_all_filt3, "/mnt/beegfs6/home3/reid/av638/ENCODE/pairwise/intersected_one_all_filt3.txt", sep="\t", append = F, quote = F, row.names = F)

intersected_one_all_re <- intersected_one_all_filt3[,c(1,2,11)]
colnames(intersected_one_all_re) <- c("CAP_a", "CAP_b", "score")
write.table(intersected_one_all_re, "/mnt/beegfs6/home3/reid/av638/ENCODE/pairwise/intersected_one_all_re.txt", sep="\t", append = F, quote = F, row.names = F)
hist(intersected_one_all_re$score, breaks = 100)

intersected_one_all_re_tab <- tidyr::pivot_wider(intersected_one_all_re, names_from = CAP_b, values_from = score)

#NA will come for those with no combinations matched
intersected_one_all_re_df <- as.data.frame(intersected_one_all_re_tab)
rownames(intersected_one_all_re_df) <- intersected_one_all_re_df$CAP_a
intersected_one_all_re_df <- intersected_one_all_re_df[,-1]
dim(intersected_one_all_re_df)

#replace NA with 0 to show the data was obsent
intersected_one_all_re_df[is.na(intersected_one_all_re_df)] <- 0

dim(intersected_one_all_re_df)

#write.table(intersected_one_all_re_df, "intersected_one_all_re_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

breaksListd <- seq(0, 1, by = 0.01)
library(pheatmap)
pheatmap_intersected_one_all_re_df <- pheatmap(as.matrix(intersected_one_all_re_df), 
                                                      na_col = "grey",
                                                      color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
                                                      clustering_distance_cols = "euclidean",
                                                      cluster_rows = T,
                                                      cluster_cols = T,
                                                      clustering_method = "ward.D")
#save as pheatmap_intersected_one_all_re_df.pdf



# summit - based






shell_summit_intersector_output <- fread("/mnt/home3/reid/av638/ENCODE/summits/shell_summit_intersector_output.txt")
num_proteins = 854
summit_mydata_p1 <- data.frame(shell_summit_intersector_output)
colnames(summit_mydata_p1) <- c("proteinA", "proteinB", "protein_Apeaks", "protein_Bpeaks", "Overlapspeaks")
summit_mydata_p1$proteinA <- gsub("_summits_sorted_merged.bed", "", summit_mydata_p1$proteinA)
summit_mydata_p1$proteinB <- gsub("_summits_sorted_merged.bed", "",summit_mydata_p1$proteinB)
summit_mydata_p1["new_expected_value"] <- round((as.numeric(summit_mydata_p1$protein_Apeaks) * as.numeric(summit_mydata_p1$protein_Bpeaks))/250000) #lambda(expected_value)=(n*m)/N
summit_mydata_p1["pois_pval"] <- dpois(summit_mydata_p1$Overlapspeaks, summit_mydata_p1$new_expected_value)
#Perform Bonferroni Correction
summit_mydata_p1["Corrected_pois_pval"] <- p.adjust(summit_mydata_p1$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
summit_mydata_p1["Corrected_pois_pval"] <- ifelse(summit_mydata_p1$Corrected_pois_pval == 0, summit_mydata_p1$Corrected_pois_pval + 1e-300, summit_mydata_p1$Corrected_pois_pval)
summit_mydata_p1["Coassociation_score"] <- ifelse((summit_mydata_p1$Overlapspeaks > summit_mydata_p1$new_expected_value), -log(summit_mydata_p1$Corrected_pois_pval,10), log(summit_mydata_p1$Corrected_pois_pval,10))
summit_mydata_p1 <- summit_mydata_p1 %>% distinct()
summit_mydata_p1_re <- summit_mydata_p1[,c(1,2,9)]
colnames(summit_mydata_p1_re) <- c("protein1", "protein2", "score")


hist(summit_mydata_p1_re$score, breaks = 100)

summit_mydata_p1_re_tab <- tidyr::pivot_wider(summit_mydata_p1_re, names_from = protein2, values_from = score)
summit_mydata_p1_re_df <- as.data.frame(summit_mydata_p1_re_tab)
rownames(summit_mydata_p1_re_df) <- summit_mydata_p1_re_df$protein1
summit_mydata_p1_re_df <- summit_mydata_p1_re_df[,-1]
head(summit_mydata_p1_re_df)

dim(summit_mydata_p1_re_df)

summit_mydata_p1_re_df_filt <- summit_mydata_p1_re_df
summit_mydata_p1_re_df_filt[is.na(summit_mydata_p1_re_df_filt)] <- 0


dim(summit_mydata_p1_re_df_filt)


# write.table(summit_mydata_p1_re_df, "summit_mydata_p1_re_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)
# write.csv(summit_mydata_p1_re_df, "summit_mydata_p1_re_df.csv")
# 
# write.table(summit_mydata_p1_re_df_filt, "summit_mydata_p1_re_df_filt.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)
# write.csv(summit_mydata_p1_re_df_filt, "summit_mydata_p1_re_df_filt.csv")

breaksListd <- seq(-300, 300, by = 0.01)
library(pheatmap)
pheatmap_summit_mydata_p1_re_df_filt <- pheatmap(as.matrix(summit_mydata_p1_re_df_filt), 
                                                 na_col = "grey",
                                                 color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
                                                 clustering_distance_cols = "correlation",
                                                 cluster_rows = T,
                                                 cluster_cols = T,
                                                 clustering_method = "ward.D")
#save as pheatmap_summit_mydata_p1_re_df_filt.pdf


# prepare hierarchical cluster
summit_mydata_p1_re_df_filt_hc = hclust(dist(t(summit_mydata_p1_re_df_filt)))
plot(summit_mydata_p1_re_df_filt_hc)



library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)


num_proteins = 329
test_p1 <- fread("/mnt/home3/reid/av638/ENCODE/summits/test.table2.txt")
test_p1 <- data.frame(test_p1)
test_p1["new_expected_value"] <- round((test_p1$TF_A_peaks * test_p1$TF_B_peaks)/250000) #lambda(expected_value)=(n*m)/N
test_p1["pois_pval"] <- dpois(test_p1$Overlapping_peaks, test_p1$new_expected_value)
#Perform Bonferroni Correction
test_p1["Corrected_pois_pval"] <- p.adjust(test_p1$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
test_p1["Corrected_pois_pval"] <- ifelse(test_p1$Corrected_pois_pval == 0, test_p1$Corrected_pois_pval + 1e-300, test_p1$Corrected_pois_pval)
test_p1["Coassociation_score"] <- ifelse((test_p1$Overlapping_peaks > test_p1$expected_value), -log(test_p1$Corrected_pois_pval,10), log(test_p1$Corrected_pois_pval,10))
test_p1 <- test_p1 %>% distinct()
test_p1_re <- test_p1[,c(4,5,15)]
colnames(test_p1_re) <- c("protein1", "protein2", "score")
test_p1_re_tab <- tidyr::pivot_wider(test_p1_re, names_from = protein2, values_from = score)
test_p1_re_df <- as.data.frame(test_p1_re_tab)
rownames(test_p1_re_df) <- test_p1_re_df$protein1
test_p1_re_df <- test_p1_re_df[,-1]
head(test_p1_re_df)

dim(test_p1_re_df)

test_p1_re_df_filt <- test_p1_re_df
test_p1_re_df_filt[is.na(test_p1_re_df_filt)] <- 0

dim(test_p1_re_df)

# write.table(test_p1_re_df, "test_p1_re_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)
# write.csv(test_p1_re_df, "test_p1_re_df.csv")
# 
# write.table(test_p1_re_df_filt, "test_p1_re_df_filt.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)
# write.csv(test_p1_re_df_filt, "test_p1_re_df_filt.csv")

breaksListd <- seq(-300, 300, by = 0.01)
library(pheatmap)
pheatmap_test_p1_re_df_filt <- pheatmap(as.matrix(test_p1_re_df_filt), 
                                        na_col = "grey",
                                        color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
                                        clustering_distance_cols = "correlation",
                                        cluster_rows = T,
                                        cluster_cols = T,
                                        clustering_method = "ward.D")
#save as pheatmap_test_p1_re_df_filt.pdf


summit_mydata_p1_hsr_re.gt100 <- fread("/mnt/home3/reid/av638/tutorial/clustering/summit_mydata_p1_hsr_re.gt100.txt")
summit_mydata_p1_hsr_re.gt100 <- data.frame(summit_mydata_p1_hsr_re.gt100)
test_p1_re <- summit_mydata_p1_hsr_re.gt100
colnames(test_p1_re) <- c("protein1", "protein2", "score")
test_p1_re_tab <- tidyr::pivot_wider(test_p1_re, names_from = protein2, values_from = score)
test_p1_re_df <- as.data.frame(test_p1_re_tab)
rownames(test_p1_re_df) <- test_p1_re_df$protein1
test_p1_re_df <- test_p1_re_df[,-1]
head(test_p1_re_df)

dim(test_p1_re_df)

test_p1_re_df_filt <- test_p1_re_df
test_p1_re_df_filt[is.na(test_p1_re_df_filt)] <- 0
# prepare hierarchical cluster
test_p1_re_df_filt_hc = hclust(dist(t(test_p1_re_df_filt)))
plot(test_p1_re_df_filt_hc)
# test_p1_re_df_filt_hc.pdf



# ---------------------
num_proteins = 852
summit_mydata_p1_hsr <- read.table("/mnt/home3/reid/av638/ENCODE/hotspotfree/shell_summit_hsr1000_50_intersector_output.txt", header = F, stringsAsFactors = F)
summit_mydata_p1_hsr <- data.frame(summit_mydata_p1_hsr)
colnames(summit_mydata_p1_hsr) <- c("proteinA", "proteinB", "protein_Apeaks", "protein_Bpeaks", "Overlapspeaks")
summit_mydata_p1_hsr$proteinA <- gsub("_summits_1000hsr50_sorted_merged.bed", "", summit_mydata_p1_hsr$proteinA)
summit_mydata_p1_hsr$proteinB <- gsub("_summits_1000hsr50_sorted_merged.bed", "",summit_mydata_p1_hsr$proteinB)
summit_mydata_p1_hsr["CAP_a"] <- unlist(lapply(strsplit(summit_mydata_p1_hsr$proteinA, "_"), function(x) x[1]))
summit_mydata_p1_hsr["CAP_b"] <- unlist(lapply(strsplit(summit_mydata_p1_hsr$proteinB, "_"), function(x) x[1]))
summit_mydata_p1_hsr["id"] <- paste0(summit_mydata_p1_hsr$CAP_a, "%" , summit_mydata_p1_hsr$CAP_b)

summit_mydata_p1_hsr_agg = aggregate(summit_mydata_p1_hsr[,c(3,4,5)], by=list(summit_mydata_p1_hsr$id), mean)
summit_mydata_p1_hsr_agg$Overlapspeaks <- round(summit_mydata_p1_hsr_agg$Overlapspeaks)
summit_mydata_p1_hsr_agg$protein_Apeaks <- round(summit_mydata_p1_hsr_agg$protein_Apeaks)
summit_mydata_p1_hsr_agg$protein_Bpeaks <- round(summit_mydata_p1_hsr_agg$protein_Bpeaks)
summit_mydata_p1_hsr_agg1 <- data.frame(cSplit(summit_mydata_p1_hsr_agg, "Group.1", "%"))
summit_mydata_p1_hsr_mg <- summit_mydata_p1_hsr_agg1[,c(4,5,1,2,3)]
summit_mydata_p1_hsr_mg["new_expected_value"] <- round((as.numeric(summit_mydata_p1_hsr_mg$protein_Apeaks) * as.numeric(summit_mydata_p1_hsr_mg$protein_Bpeaks))/250000) #lambda(expected_value)=(n*m)/N
summit_mydata_p1_hsr_mg["pois_pval"] <- dpois(summit_mydata_p1_hsr_mg$Overlapspeaks, summit_mydata_p1_hsr_mg$new_expected_value)
#Perform Bonferroni Correction
summit_mydata_p1_hsr_mg["Corrected_pois_pval"] <- p.adjust(summit_mydata_p1_hsr_mg$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
summit_mydata_p1_hsr_mg["Corrected_pois_pval"] <- ifelse(summit_mydata_p1_hsr_mg$Corrected_pois_pval == 0, summit_mydata_p1_hsr_mg$Corrected_pois_pval + 1e-300, summit_mydata_p1_hsr_mg$Corrected_pois_pval)
summit_mydata_p1_hsr_mg["Coassociation_score"] <- ifelse((summit_mydata_p1_hsr_mg$Overlapspeaks > summit_mydata_p1_hsr_mg$new_expected_value), -log(summit_mydata_p1_hsr_mg$Corrected_pois_pval,10), log(summit_mydata_p1_hsr_mg$Corrected_pois_pval,10))
summit_mydata_p1_hsr_mg <- summit_mydata_p1_hsr_mg %>% distinct()

summit_mydata_p1_hsr_mg["proteinA"] <- summit_mydata_p1_hsr_mg$Group.1_1
summit_mydata_p1_hsr_mg["proteinB"] <- summit_mydata_p1_hsr_mg$Group.1_2
summit_mydata_p1_hsr_mg <- summit_mydata_p1_hsr_mg[,c(10,11,3:9)]
#filter low count peaks
#sideA
summit_mydata_p1_hsr_mg_filt <- summit_mydata_p1_hsr_mg[which(summit_mydata_p1_hsr_mg$protein_Apeaks > 100),]
#sideB
summit_mydata_p1_hsr_mg_filt <- summit_mydata_p1_hsr_mg_filt[which(summit_mydata_p1_hsr_mg_filt$protein_Bpeaks > 100),]

#remove self-hits same protein on both sides
summit_mydata_p1_hsr_mg_filt <- summit_mydata_p1_hsr_mg_filt[which(summit_mydata_p1_hsr_mg_filt$proteinA != summit_mydata_p1_hsr_mg_filt$proteinB),]

#remove 0s coassociation score
summit_mydata_p1_hsr_mg_filt2 <- summit_mydata_p1_hsr_mg_filt[which(summit_mydata_p1_hsr_mg_filt$Coassociation_score != 0),]

#remove +/- 1.4 coassociation
summit_mydata_p1_hsr_mg_filt3 <- summit_mydata_p1_hsr_mg_filt2[which(summit_mydata_p1_hsr_mg_filt2$Coassociation_score > 1.4 | summit_mydata_p1_hsr_mg_filt2$Coassociation_score < -1.4),]
#filter non-significant 0.05
# for (i in unique(summit_mydata_p1_hsr_mg$proteinA)){
#   for (j in unique(summit_mydata_p1_hsr_mg$proteinB)){
#     if (i != j){
#       print(paste0(i ,j))
#     }
#   }
# }

# write.table(summit_mydata_p1_hsr_mg_filt3, "summit_mydata_p1_hsr_mg_filt3.txt", sep="\t", append = F, quote = F, row.names = F)

summit_mydata_p1_hsr_mg_re <- summit_mydata_p1_hsr_mg_filt3[,c(1,2,9)]
colnames(summit_mydata_p1_hsr_mg_re) <- c("protein1", "protein2", "score")
hist(summit_mydata_p1_hsr_mg_re$score, breaks = 100)

write.table(summit_mydata_p1_hsr_mg_re, "/mnt/home3/reid/av638/tutorial/clustering/summit_mydata_p1_hsr_mg_re.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)


summit_mydata_p1_hsr_mg_re.gt100 <- summit_mydata_p1_hsr_mg_re[which(summit_mydata_p1_hsr_mg_re$score >= 100),]
write.table(summit_mydata_p1_hsr_mg_re.gt100, "/mnt/home3/reid/av638/tutorial/clustering/summit_mydata_p1_hsr_mg_re.gt100.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)

summit_mydata_p1_hsr_mg_re.gt10 <- summit_mydata_p1_hsr_mg_re[which(summit_mydata_p1_hsr_mg_re$score >= 10),]
write.table(summit_mydata_p1_hsr_mg_re.gt10, "/mnt/home3/reid/av638/tutorial/clustering/summit_mydata_p1_hsr_mg_re.gt10.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)


# #remove some samples
# samples_to_remove <- read.table("samples_to_remove.txt")
# 
# list_samples_to_keep <- data.frame("ids"=setdiff(summit_mydata_p1_hsr_mg_re$protein1, samples_to_remove$V1))
# summit_mydata_p1_hsr_mg_retain1 <- merge(summit_mydata_p1_hsr_mg_re, list_samples_to_keep, by.x ="protein1", by.y="ids")
# summit_mydata_p1_hsr_mg_retain2 <- merge(summit_mydata_p1_hsr_mg_retain1, list_samples_to_keep, by.x ="protein2", by.y="ids")
# 
# write.table(summit_mydata_p1_hsr_mg_retain2, "summit_mydata_p1_hsr_mg_retained.txt", sep="\t", append = F, quote = F, row.names = F)

summit_mydata_p1_hsr_mg_re_tab <- tidyr::pivot_wider(summit_mydata_p1_hsr_mg_re, names_from = protein2, values_from = score)

#NA will come for those with no combinations matched
summit_mydata_p1_hsr_mg_re_df <- as.data.frame(summit_mydata_p1_hsr_mg_re_tab)
rownames(summit_mydata_p1_hsr_mg_re_df) <- summit_mydata_p1_hsr_mg_re_df$protein1
summit_mydata_p1_hsr_mg_re_df <- summit_mydata_p1_hsr_mg_re_df[,-1]
dim(summit_mydata_p1_hsr_mg_re_df)

#replace NA with 0 to show the data was obsent
summit_mydata_p1_hsr_mg_re_df[is.na(summit_mydata_p1_hsr_mg_re_df)] <- 0

dim(summit_mydata_p1_hsr_mg_re_df)

write.table(summit_mydata_p1_hsr_mg_re_df, "/mnt/home3/reid/av638/tutorial/clustering/summit_mydata_p1_hsr_mg_re_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

breaksListd <- seq(0, 1, by = 0.01)
library(pheatmap)
pheatmap_summit_mydata_p1_hsr_mg_re_df <- pheatmap(summit_mydata_p1_hsr_mg_re_df, 
                                                   na_col = "grey",
                                                   color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
                                                   clustering_distance_cols = "euclidean",
                                                   cluster_rows = T,
                                                   cluster_cols = T,
                                                   clustering_method = "ward.D")


# prepare hierarchical cluster
summit_mydata_p1_hsr_mg_re_df_hc = hclust(dist(t(summit_mydata_p1_hsr_mg_re_df)))
plot(summit_mydata_p1_hsr_mg_re_df_hc)
# summit_mydata_p1_hsr_mg_re_df_hc.pdf

# separate groups
summit_mydata_p1_hsr_mg_re_df_cutree <- data.frame(groups=cutree(summit_mydata_p1_hsr_mg_re_df_hc, k = 50))
summit_mydata_p1_hsr_mg_re_df_cutree["caps"] <- rownames(summit_mydata_p1_hsr_mg_re_df_cutree)
count_summit_mydata_p1_hsr_mg_re_df_cutree <- count(summit_mydata_p1_hsr_mg_re_df_cutree$groups)
colnames(count_summit_mydata_p1_hsr_mg_re_df_cutree) <- c("complex", "num_of_caps")
plot(count_summit_mydata_p1_hsr_mg_re_df_cutree$complex, count_summit_mydata_p1_hsr_mg_re_df_cutree$num_of_caps)

summit_mydata_p1_hsr_mg_re_df_cutree_merge <- ddply(summit_mydata_p1_hsr_mg_re_df_cutree, .(groups), summarize, caps = toString(caps))
# summit_mydata_p1_hsr_mg_re_df_cutree_merge[summit_mydata_p1_hsr_mg_re_df_cutree_merge$caps %like% "CTCF", ]
summit_mydata_p1_hsr_mg_re_df_cutree_merge[grepl("\\CTCF\\b", summit_mydata_p1_hsr_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]
summit_mydata_p1_hsr_mg_re_df_cutree_merge[grepl("\\SUZ12\\b", summit_mydata_p1_hsr_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]
summit_mydata_p1_hsr_mg_re_df_cutree_merge[grepl("\\NFYA\\b", summit_mydata_p1_hsr_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]

# use graph based clustering to identify complexes
# Cluster data with MCL
mcxload -abc summit_mydata_p1_hsr_re.gt100.txt --stream-mirror -write-tab data.tab -o data.mci 
# mcxload -abc summit_mydata_p1_hsr_mg_re.gt100.txt --stream-mirror -write-tab data.tab -o data.mci 
# mcxload -abc summit_mydata_p1_hsr_mg_re.gt10.txt --stream-mirror -write-tab data.tab -o data.mci 


mcl data.mci -I 1.4
mcl data.mci -I 2
mcl data.mci -I 4

mcxdump -icl out.data.mci.I14 -tabr data.tab -o dump.data.mci.I14
mcxdump -icl out.data.mci.I20 -tabr data.tab -o dump.data.mci.I20
mcxdump -icl out.data.mci.I40 -tabr data.tab -o dump.data.mci.I40

# Look at cluster sizes
cat dump.data.mci.I14  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I20  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I40  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'



# Obtain known protein complex data from  http://mips.helmholtz-muenchen.de/corum/#download
  
# Filter known complexes to only contain proteins from the ENCODE dataset
python complex_count_filter.py humanComplexes.txt summit_mydata_p1_hsr_re.txt > humanComplexes_filtered.txt

# Look at accuracy of predicted complexes
python mod_benchmark_complexes.py dump.data.mci.I14 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'

# Which complexes do we predict perfectly?
python mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'

# this script plots the senstivity and specificity to show improvements as inflation parameter increases
complex_analysis.R

library(readxl)
# http://human.med.utoronto.ca/php/data_download.php
# cite: https://www.cell.com/fulltext/S0092-8674(12)01006-9
human.med.utoronto.ca_TableS3 <- readxl::read_excel("/mnt/home3/reid/av638/tutorial/clustering/human.med.utoronto.ca_TableS3.xls", sheet = 1, col_names = FALSE, skip=1)
human.med.utoronto.ca_TableS3 <- data.frame(human.med.utoronto.ca_TableS3)
human.med.utoronto.ca_complex <- data.frame(complex=human.med.utoronto.ca_TableS3$...4)
human.med.utoronto.ca_complex <- data.frame(complex=human.med.utoronto.ca_complex[-1,])

human.med.utoronto.ca_complex_list <- lapply(strsplit(human.med.utoronto.ca_complex$complex, ","), function(x) paste0(gsub("_HUMAN","", x), collapse = ","))
human.med.utoronto.ca_complex_df <- do.call(rbind.data.frame, human.med.utoronto.ca_complex_list)
colnames(human.med.utoronto.ca_complex_df) <- "CAPs"
human.med.utoronto.ca_complex_df["complex_id"] <- paste0("complex_",rownames(human.med.utoronto.ca_complex_df))
human.med.utoronto.ca_complex_df <- human.med.utoronto.ca_complex_df[,c(2,1)]
write.table(human.med.utoronto.ca_complex_df, "/mnt/home3/reid/av638/tutorial/clustering/human.med.utoronto.ca_complex_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

# python benchmarkagain_complex.py dump.data.mci.I40 human.med.utoronto.ca_complex_df.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# python benchmarkagain_complex.py dump.data.mci.I40 human.med.utoronto.ca_complex_df.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'


# http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt
# cite:  https://www.embopress.org/doi/full/10.15252/msb.20167490

humap2_complexes_20200809 <- fread("/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809.txt")
humap2_complexes_20200809 <- data.frame(humap2_complexes_20200809)
humap2_complexes_20200809["genes"] <- unlist(lapply(strsplit(humap2_complexes_20200809$genenames, " "), function(x) paste0(x, collapse = ",")))
humap2_complexes_20200809_re <- humap2_complexes_20200809[,c(1,5)]
write.table(humap2_complexes_20200809_re, "/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809_re.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

# python benchmarkagain_complex.py dump.data.mci.I40 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# python benchmarkagain_complex.py dump.data.mci.I40 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'




# ---------------------
num_proteins = 852
summit_mydata_p1_cbh <- read.table("/mnt/beegfs6/home3/reid/av638/ENCODE/summit_pairwise/shell_summit_qualblackhot1000_50_intersector_output.txt", header = F, stringsAsFactors = F)
summit_mydata_p1_cbh <- data.frame(summit_mydata_p1_cbh)
colnames(summit_mydata_p1_cbh) <- c("proteinA", "proteinB", "protein_Apeaks", "protein_Bpeaks", "Overlapspeaks")
summit_mydata_p1_cbh$proteinA <- gsub("_summits_1000hsr50_sorted_merged.bed", "", summit_mydata_p1_cbh$proteinA)
summit_mydata_p1_cbh$proteinB <- gsub("_summits_1000hsr50_sorted_merged.bed", "",summit_mydata_p1_cbh$proteinB)
summit_mydata_p1_cbh["CAP_a"] <- unlist(lapply(strsplit(summit_mydata_p1_cbh$proteinA, "_"), function(x) x[1]))
summit_mydata_p1_cbh["CAP_b"] <- unlist(lapply(strsplit(summit_mydata_p1_cbh$proteinB, "_"), function(x) x[1]))
summit_mydata_p1_cbh["id"] <- paste0(summit_mydata_p1_cbh$CAP_a, "%" , summit_mydata_p1_cbh$CAP_b)

summit_mydata_p1_cbh_agg = aggregate(summit_mydata_p1_cbh[,c(3,4,5)], by=list(summit_mydata_p1_cbh$id), mean)
summit_mydata_p1_cbh_agg$Overlapspeaks <- round(summit_mydata_p1_cbh_agg$Overlapspeaks)
summit_mydata_p1_cbh_agg$protein_Apeaks <- round(summit_mydata_p1_cbh_agg$protein_Apeaks)
summit_mydata_p1_cbh_agg$protein_Bpeaks <- round(summit_mydata_p1_cbh_agg$protein_Bpeaks)
summit_mydata_p1_cbh_agg1 <- data.frame(cSplit(summit_mydata_p1_cbh_agg, "Group.1", "%"))
summit_mydata_p1_cbh_mg <- summit_mydata_p1_cbh_agg1[,c(4,5,1,2,3)]
summit_mydata_p1_cbh_mg["new_expected_value"] <- round((as.numeric(summit_mydata_p1_cbh_mg$protein_Apeaks) * as.numeric(summit_mydata_p1_cbh_mg$protein_Bpeaks))/250000) #lambda(expected_value)=(n*m)/N
summit_mydata_p1_cbh_mg["pois_pval"] <- dpois(summit_mydata_p1_cbh_mg$Overlapspeaks, summit_mydata_p1_cbh_mg$new_expected_value)
#Perform Bonferroni Correction
summit_mydata_p1_cbh_mg["Corrected_pois_pval"] <- p.adjust(summit_mydata_p1_cbh_mg$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
summit_mydata_p1_cbh_mg["Corrected_pois_pval"] <- ifelse(summit_mydata_p1_cbh_mg$Corrected_pois_pval == 0, summit_mydata_p1_cbh_mg$Corrected_pois_pval + 1e-300, summit_mydata_p1_cbh_mg$Corrected_pois_pval)
summit_mydata_p1_cbh_mg["Coassociation_score"] <- ifelse((summit_mydata_p1_cbh_mg$Overlapspeaks > summit_mydata_p1_cbh_mg$new_expected_value), -log(summit_mydata_p1_cbh_mg$Corrected_pois_pval,10), log(summit_mydata_p1_cbh_mg$Corrected_pois_pval,10))
summit_mydata_p1_cbh_mg <- summit_mydata_p1_cbh_mg %>% distinct()

summit_mydata_p1_cbh_mg["proteinA"] <- summit_mydata_p1_cbh_mg$Group.1_1
summit_mydata_p1_cbh_mg["proteinB"] <- summit_mydata_p1_cbh_mg$Group.1_2
summit_mydata_p1_cbh_mg <- summit_mydata_p1_cbh_mg[,c(10,11,3:9)]
#filter low count peaks
#sideA
summit_mydata_p1_cbh_mg_filt <- summit_mydata_p1_cbh_mg[which(summit_mydata_p1_cbh_mg$protein_Apeaks > 100),]
#sideB
summit_mydata_p1_cbh_mg_filt <- summit_mydata_p1_cbh_mg_filt[which(summit_mydata_p1_cbh_mg_filt$protein_Bpeaks > 100),]

#remove self-hits same protein on both sides
summit_mydata_p1_cbh_mg_filt <- summit_mydata_p1_cbh_mg_filt[which(summit_mydata_p1_cbh_mg_filt$proteinA != summit_mydata_p1_cbh_mg_filt$proteinB),]

#remove 0s coassociation score
summit_mydata_p1_cbh_mg_filt2 <- summit_mydata_p1_cbh_mg_filt[which(summit_mydata_p1_cbh_mg_filt$Coassociation_score != 0),]

#remove +/- 1.4 coassociation
summit_mydata_p1_cbh_mg_filt3 <- summit_mydata_p1_cbh_mg_filt2[which(summit_mydata_p1_cbh_mg_filt2$Coassociation_score > 1.4 | summit_mydata_p1_cbh_mg_filt2$Coassociation_score < -1.4),]
#filter non-significant 0.05
# for (i in unique(summit_mydata_p1_cbh_mg$proteinA)){
#   for (j in unique(summit_mydata_p1_cbh_mg$proteinB)){
#     if (i != j){
#       print(paste0(i ,j))
#     }
#   }
# }

# write.table(summit_mydata_p1_cbh_mg_filt3, "summit_mydata_p1_cbh_mg_filt3.txt", sep="\t", append = F, quote = F, row.names = F)

summit_mydata_p1_cbh_mg_re <- summit_mydata_p1_cbh_mg_filt3[,c(1,2,9)]
colnames(summit_mydata_p1_cbh_mg_re) <- c("protein1", "protein2", "score")
hist(summit_mydata_p1_cbh_mg_re$score, breaks = 100)

write.table(summit_mydata_p1_cbh_mg_re, "/mnt/home3/reid/av638/tutorial/clustering/summit_mydata_p1_cbh_mg_re.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)


summit_mydata_p1_cbh_mg_re.gt100 <- summit_mydata_p1_cbh_mg_re[which(summit_mydata_p1_cbh_mg_re$score >= 100),]
write.table(summit_mydata_p1_cbh_mg_re.gt100, "/mnt/home3/reid/av638/tutorial/clustering/summit_mydata_p1_cbh_mg_re.gt100.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)

summit_mydata_p1_cbh_mg_re.gt10 <- summit_mydata_p1_cbh_mg_re[which(summit_mydata_p1_cbh_mg_re$score >= 10),]
write.table(summit_mydata_p1_cbh_mg_re.gt10, "/mnt/home3/reid/av638/tutorial/clustering/summit_mydata_p1_cbh_mg_re.gt10.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)


# #remove some samples
# samples_to_remove <- read.table("samples_to_remove.txt")
# 
# list_samples_to_keep <- data.frame("ids"=setdiff(summit_mydata_p1_cbh_mg_re$protein1, samples_to_remove$V1))
# summit_mydata_p1_cbh_mg_retain1 <- merge(summit_mydata_p1_cbh_mg_re, list_samples_to_keep, by.x ="protein1", by.y="ids")
# summit_mydata_p1_cbh_mg_retain2 <- merge(summit_mydata_p1_cbh_mg_retain1, list_samples_to_keep, by.x ="protein2", by.y="ids")
# 
# write.table(summit_mydata_p1_cbh_mg_retain2, "summit_mydata_p1_cbh_mg_retained.txt", sep="\t", append = F, quote = F, row.names = F)

summit_mydata_p1_cbh_mg_re_tab <- tidyr::pivot_wider(summit_mydata_p1_cbh_mg_re, names_from = protein2, values_from = score)

#NA will come for those with no combinations matched
summit_mydata_p1_cbh_mg_re_df <- as.data.frame(summit_mydata_p1_cbh_mg_re_tab)
rownames(summit_mydata_p1_cbh_mg_re_df) <- summit_mydata_p1_cbh_mg_re_df$protein1
summit_mydata_p1_cbh_mg_re_df <- summit_mydata_p1_cbh_mg_re_df[,-1]
dim(summit_mydata_p1_cbh_mg_re_df)

#replace NA with 0 to show the data was obsent
summit_mydata_p1_cbh_mg_re_df[is.na(summit_mydata_p1_cbh_mg_re_df)] <- 0

dim(summit_mydata_p1_cbh_mg_re_df)

write.table(summit_mydata_p1_cbh_mg_re_df, "/mnt/home3/reid/av638/tutorial/clustering/summit_mydata_p1_cbh_mg_re_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

breaksListd <- seq(0, 1, by = 0.01)
library(pheatmap)
pheatmap_summit_mydata_p1_cbh_mg_re_df <- pheatmap(summit_mydata_p1_cbh_mg_re_df, 
                                                   na_col = "grey",
                                                   color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
                                                   clustering_distance_cols = "euclidean",
                                                   cluster_rows = T,
                                                   cluster_cols = T,
                                                   clustering_method = "ward.D")


# prepare hierarchical cluster
summit_mydata_p1_cbh_mg_re_df_hc = hclust(dist(t(summit_mydata_p1_cbh_mg_re_df)))
plot(summit_mydata_p1_cbh_mg_re_df_hc)
# summit_mydata_p1_cbh_mg_re_df_hc.pdf

# separate groups
summit_mydata_p1_cbh_mg_re_df_cutree <- data.frame(groups=cutree(summit_mydata_p1_cbh_mg_re_df_hc, k = 50))
summit_mydata_p1_cbh_mg_re_df_cutree["caps"] <- rownames(summit_mydata_p1_cbh_mg_re_df_cutree)
count_summit_mydata_p1_cbh_mg_re_df_cutree <- count(summit_mydata_p1_cbh_mg_re_df_cutree$groups)
colnames(count_summit_mydata_p1_cbh_mg_re_df_cutree) <- c("complex", "num_of_caps")
plot(count_summit_mydata_p1_cbh_mg_re_df_cutree$complex, count_summit_mydata_p1_cbh_mg_re_df_cutree$num_of_caps)

summit_mydata_p1_cbh_mg_re_df_cutree_merge <- ddply(summit_mydata_p1_cbh_mg_re_df_cutree, .(groups), summarize, caps = toString(caps))
# summit_mydata_p1_cbh_mg_re_df_cutree_merge[summit_mydata_p1_cbh_mg_re_df_cutree_merge$caps %like% "CTCF", ]
summit_mydata_p1_cbh_mg_re_df_cutree_merge[grepl("\\CTCF\\b", summit_mydata_p1_cbh_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]
summit_mydata_p1_cbh_mg_re_df_cutree_merge[grepl("\\SUZ12\\b", summit_mydata_p1_cbh_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]
summit_mydata_p1_cbh_mg_re_df_cutree_merge[grepl("\\NFYA\\b", summit_mydata_p1_cbh_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]

# use graph based clustering to identify complexes
# Cluster data with MCL
mcxload -abc summit_mydata_p1_cbh_re.gt100.txt --stream-mirror -write-tab data.tab -o data.mci 
# mcxload -abc summit_mydata_p1_cbh_mg_re.gt100.txt --stream-mirror -write-tab data.tab -o data.mci 
# mcxload -abc summit_mydata_p1_cbh_mg_re.gt10.txt --stream-mirror -write-tab data.tab -o data.mci 


mcl data.mci -I 1.4
mcl data.mci -I 2
mcl data.mci -I 4

mcxdump -icl out.data.mci.I14 -tabr data.tab -o dump.data.mci.I14
mcxdump -icl out.data.mci.I20 -tabr data.tab -o dump.data.mci.I20
mcxdump -icl out.data.mci.I40 -tabr data.tab -o dump.data.mci.I40

# Look at cluster sizes
cat dump.data.mci.I14  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I20  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I40  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'



# Obtain known protein complex data from  http://mips.helmholtz-muenchen.de/corum/#download

# Filter known complexes to only contain proteins from the ENCODE dataset
python complex_count_filter.py humanComplexes.txt summit_mydata_p1_cbh_re.txt > humanComplexes_filtered.txt

# Look at accuracy of predicted complexes
python mod_benchmark_complexes.py dump.data.mci.I14 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'

# Which complexes do we predict perfectly?
python mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'

# this script plots the senstivity and specificity to show improvements as inflation parameter increases
complex_analysis.R

library(readxl)
# http://human.med.utoronto.ca/php/data_download.php
# cite: https://www.cell.com/fulltext/S0092-8674(12)01006-9
human.med.utoronto.ca_TableS3 <- readxl::read_excel("/mnt/home3/reid/av638/tutorial/clustering/human.med.utoronto.ca_TableS3.xls", sheet = 1, col_names = FALSE, skip=1)
human.med.utoronto.ca_TableS3 <- data.frame(human.med.utoronto.ca_TableS3)
human.med.utoronto.ca_complex <- data.frame(complex=human.med.utoronto.ca_TableS3$...4)
human.med.utoronto.ca_complex <- data.frame(complex=human.med.utoronto.ca_complex[-1,])

human.med.utoronto.ca_complex_list <- lapply(strsplit(human.med.utoronto.ca_complex$complex, ","), function(x) paste0(gsub("_HUMAN","", x), collapse = ","))
human.med.utoronto.ca_complex_df <- do.call(rbind.data.frame, human.med.utoronto.ca_complex_list)
colnames(human.med.utoronto.ca_complex_df) <- "CAPs"
human.med.utoronto.ca_complex_df["complex_id"] <- paste0("complex_",rownames(human.med.utoronto.ca_complex_df))
human.med.utoronto.ca_complex_df <- human.med.utoronto.ca_complex_df[,c(2,1)]
write.table(human.med.utoronto.ca_complex_df, "/mnt/home3/reid/av638/tutorial/clustering/human.med.utoronto.ca_complex_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

# python benchmarkagain_complex.py dump.data.mci.I40 human.med.utoronto.ca_complex_df.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# python benchmarkagain_complex.py dump.data.mci.I40 human.med.utoronto.ca_complex_df.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'


# http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt
# cite:  https://www.embopress.org/doi/full/10.15252/msb.20167490

humap2_complexes_20200809 <- fread("/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809.txt")
humap2_complexes_20200809 <- data.frame(humap2_complexes_20200809)
humap2_complexes_20200809["genes"] <- unlist(lapply(strsplit(humap2_complexes_20200809$genenames, " "), function(x) paste0(x, collapse = ",")))
humap2_complexes_20200809_re <- humap2_complexes_20200809[,c(1,5)]
write.table(humap2_complexes_20200809_re, "/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809_re.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

# python benchmarkagain_complex.py dump.data.mci.I40 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# python benchmarkagain_complex.py dump.data.mci.I40 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'


# ---------------------  main ------------------------- #
num_proteins = 514
summit_mydata_p1_cbh <- read.table("/mnt/beegfs6/home3/reid/av638/ENCODE/summit_pairwise/shell_summit_qualblackhot1000_50_intersector_output_all.txt", header = F, stringsAsFactors = F)
summit_mydata_p1_cbh <- data.frame(summit_mydata_p1_cbh)
colnames(summit_mydata_p1_cbh) <- c("proteinA", "proteinB", "protein_Apeaks", "protein_Bpeaks", "Overlapspeaks")


summit_mydata_p1_cbh["CAP_a"] <- summit_mydata_p1_cbh$proteinA
summit_mydata_p1_cbh["CAP_b"] <- summit_mydata_p1_cbh$proteinB
summit_mydata_p1_cbh_mg <- summit_mydata_p1_cbh[,c(6,7,3:5)]
summit_mydata_p1_cbh_mg["new_expected_value"] <- round((as.numeric(summit_mydata_p1_cbh_mg$protein_Apeaks) * as.numeric(summit_mydata_p1_cbh_mg$protein_Bpeaks))/250000) #lambda(expected_value)=(n*m)/N
summit_mydata_p1_cbh_mg["pois_pval"] <- dpois(summit_mydata_p1_cbh_mg$Overlapspeaks, summit_mydata_p1_cbh_mg$new_expected_value)
#Perform Bonferroni Correction
summit_mydata_p1_cbh_mg["Corrected_pois_pval"] <- p.adjust(summit_mydata_p1_cbh_mg$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
summit_mydata_p1_cbh_mg["Corrected_pois_pval"] <- ifelse(summit_mydata_p1_cbh_mg$Corrected_pois_pval == 0, summit_mydata_p1_cbh_mg$Corrected_pois_pval + 1e-300, summit_mydata_p1_cbh_mg$Corrected_pois_pval)
summit_mydata_p1_cbh_mg["Coassociation_score"] <- ifelse((summit_mydata_p1_cbh_mg$Overlapspeaks > summit_mydata_p1_cbh_mg$new_expected_value), -log(summit_mydata_p1_cbh_mg$Corrected_pois_pval,10), log(summit_mydata_p1_cbh_mg$Corrected_pois_pval,10))
summit_mydata_p1_cbh_mg <- summit_mydata_p1_cbh_mg %>% distinct()

summit_mydata_p1_cbh_mg["proteinA"] <- summit_mydata_p1_cbh_mg$CAP_a
summit_mydata_p1_cbh_mg["proteinB"] <- summit_mydata_p1_cbh_mg$CAP_b
summit_mydata_p1_cbh_mg <- summit_mydata_p1_cbh_mg[,c(10,11,3:9)]
#filter low count peaks
#sideA
summit_mydata_p1_cbh_mg_filt <- summit_mydata_p1_cbh_mg[which(summit_mydata_p1_cbh_mg$protein_Apeaks > 30),]
#sideB
summit_mydata_p1_cbh_mg_filt <- summit_mydata_p1_cbh_mg_filt[which(summit_mydata_p1_cbh_mg_filt$protein_Bpeaks > 30),]

#remove self-hits same protein on both sides
summit_mydata_p1_cbh_mg_filt <- summit_mydata_p1_cbh_mg_filt[which(summit_mydata_p1_cbh_mg_filt$proteinA != summit_mydata_p1_cbh_mg_filt$proteinB),]

#remove 0s coassociation score
summit_mydata_p1_cbh_mg_filt2 <- summit_mydata_p1_cbh_mg_filt[which(summit_mydata_p1_cbh_mg_filt$Coassociation_score != 0),]

#remove +/- 1.4 coassociation
summit_mydata_p1_cbh_mg_filt3 <- summit_mydata_p1_cbh_mg_filt2[which(summit_mydata_p1_cbh_mg_filt2$Coassociation_score > 1.4 | summit_mydata_p1_cbh_mg_filt2$Coassociation_score < -1.4),]
#filter non-significant 0.05
# for (i in unique(summit_mydata_p1_cbh_mg$proteinA)){
#   for (j in unique(summit_mydata_p1_cbh_mg$proteinB)){
#     if (i != j){
#       print(paste0(i ,j))
#     }
#   }
# }

# write.table(summit_mydata_p1_cbh_mg_filt3, "summit_mydata_p1_cbh_mg_filt3.txt", sep="\t", append = F, quote = F, row.names = F)

summit_mydata_p1_cbh_mg_re <- summit_mydata_p1_cbh_mg_filt3[,c(1,2,9)]
colnames(summit_mydata_p1_cbh_mg_re) <- c("protein1", "protein2", "score")
hist(summit_mydata_p1_cbh_mg_re$score, breaks = 100)

write.table(summit_mydata_p1_cbh_mg_re, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/summit_mydata_p1_cbh_mg_re.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)


summit_mydata_p1_cbh_mg_re.gt100 <- summit_mydata_p1_cbh_mg_re[which(summit_mydata_p1_cbh_mg_re$score >= 100),]
write.table(summit_mydata_p1_cbh_mg_re.gt100, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/summit_mydata_p1_cbh_mg_re.gt100.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)

summit_mydata_p1_cbh_mg_re.gt10 <- summit_mydata_p1_cbh_mg_re[which(summit_mydata_p1_cbh_mg_re$score >= 10),]
write.table(summit_mydata_p1_cbh_mg_re.gt10, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/summit_mydata_p1_cbh_mg_re.gt10.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)


summit_mydata_p1_cbh_mg_re_tab <- tidyr::pivot_wider(summit_mydata_p1_cbh_mg_re, names_from = protein2, values_from = score)

#NA will come for those with no combinations matched
summit_mydata_p1_cbh_mg_re_df <- as.data.frame(summit_mydata_p1_cbh_mg_re_tab)
rownames(summit_mydata_p1_cbh_mg_re_df) <- summit_mydata_p1_cbh_mg_re_df$protein1
summit_mydata_p1_cbh_mg_re_df <- summit_mydata_p1_cbh_mg_re_df[,-1]
dim(summit_mydata_p1_cbh_mg_re_df)

#replace NA with 0 to show the data was obsent
summit_mydata_p1_cbh_mg_re_df[is.na(summit_mydata_p1_cbh_mg_re_df)] <- 0

dim(summit_mydata_p1_cbh_mg_re_df)

write.table(summit_mydata_p1_cbh_mg_re_df, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/summit_mydata_p1_cbh_mg_re_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

breaksListd <- seq(0, 1, by = 0.01)
library(pheatmap)
pheatmap_summit_mydata_p1_cbh_mg_re_df <- pheatmap(summit_mydata_p1_cbh_mg_re_df, 
                                                   na_col = "grey",
                                                   color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
                                                   clustering_distance_cols = "euclidean",
                                                   cluster_rows = T,
                                                   cluster_cols = T,
                                                   clustering_method = "ward.D")


# prepare hierarchical cluster
summit_mydata_p1_cbh_mg_re_df_hc = hclust(dist(t(summit_mydata_p1_cbh_mg_re_df)))
plot(summit_mydata_p1_cbh_mg_re_df_hc)
# summit_mydata_p1_cbh_mg_re_df_hc.pdf

# separate groups
summit_mydata_p1_cbh_mg_re_df_cutree <- data.frame(groups=cutree(summit_mydata_p1_cbh_mg_re_df_hc, k = 50))
summit_mydata_p1_cbh_mg_re_df_cutree["caps"] <- rownames(summit_mydata_p1_cbh_mg_re_df_cutree)
count_summit_mydata_p1_cbh_mg_re_df_cutree <- count(summit_mydata_p1_cbh_mg_re_df_cutree$groups)
colnames(count_summit_mydata_p1_cbh_mg_re_df_cutree) <- c("complex", "num_of_caps")
plot(count_summit_mydata_p1_cbh_mg_re_df_cutree$complex, count_summit_mydata_p1_cbh_mg_re_df_cutree$num_of_caps)

summit_mydata_p1_cbh_mg_re_df_cutree_merge <- ddply(summit_mydata_p1_cbh_mg_re_df_cutree, .(groups), summarize, caps = toString(caps))
# summit_mydata_p1_cbh_mg_re_df_cutree_merge[summit_mydata_p1_cbh_mg_re_df_cutree_merge$caps %like% "CTCF", ]
summit_mydata_p1_cbh_mg_re_df_cutree_merge[grepl("\\CTCF\\b", summit_mydata_p1_cbh_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]
summit_mydata_p1_cbh_mg_re_df_cutree_merge[grepl("\\SUZ12\\b", summit_mydata_p1_cbh_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]
summit_mydata_p1_cbh_mg_re_df_cutree_merge[grepl("\\NFYA\\b", summit_mydata_p1_cbh_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]

# use graph based clustering to identify complexes
# Cluster data with MCL
mcxload -abc summit_mydata_p1_cbh_mg_re.gt100.txt --stream-mirror -write-tab data.tab -o data.mci 
# mcxload -abc summit_mydata_p1_cbh_mg_re.gt10.txt --stream-mirror -write-tab data.tab -o data.mci 


mcl data.mci -I 1.4
mcl data.mci -I 2
mcl data.mci -I 4

mcxdump -icl out.data.mci.I14 -tabr data.tab -o dump.data.mci.I14
mcxdump -icl out.data.mci.I20 -tabr data.tab -o dump.data.mci.I20
mcxdump -icl out.data.mci.I40 -tabr data.tab -o dump.data.mci.I40

# Look at cluster sizes
cat dump.data.mci.I14  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I20  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I40  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'



# Obtain known protein complex data from  http://mips.helmholtz-muenchen.de/corum/#download

# Filter known complexes to only contain proteins from the ENCODE dataset
python3 complex_count_filter.py humanComplexes.txt summit_mydata_p1_cbh_mg_re.txt > humanComplexes_filtered.txt

# Look at accuracy of predicted complexes
python3 mod_benchmark_complexes.py dump.data.mci.I14 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
python3 mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'

# Which complexes do we predict perfectly?
python3 mod_benchmark_complexes.py dump.data.mci.I20 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'
python3 mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'

# this script plots the senstivity and specificity to show improvements as inflation parameter increases
complex_analysis.R

library(readxl)
# http://human.med.utoronto.ca/php/data_download.php
# cite: https://www.cell.com/fulltext/S0092-8674(12)01006-9
human.med.utoronto.ca_TableS3 <- readxl::read_excel("/mnt/home3/reid/av638/tutorial/clustering/human.med.utoronto.ca_TableS3.xls", sheet = 1, col_names = FALSE, skip=1)
human.med.utoronto.ca_TableS3 <- data.frame(human.med.utoronto.ca_TableS3)
human.med.utoronto.ca_complex <- data.frame(complex=human.med.utoronto.ca_TableS3$...4)
human.med.utoronto.ca_complex <- data.frame(complex=human.med.utoronto.ca_complex[-1,])

human.med.utoronto.ca_complex_list <- lapply(strsplit(human.med.utoronto.ca_complex$complex, ","), function(x) paste0(gsub("_HUMAN","", x), collapse = ","))
human.med.utoronto.ca_complex_df <- do.call(rbind.data.frame, human.med.utoronto.ca_complex_list)
colnames(human.med.utoronto.ca_complex_df) <- "CAPs"
human.med.utoronto.ca_complex_df["complex_id"] <- paste0("complex_",rownames(human.med.utoronto.ca_complex_df))
human.med.utoronto.ca_complex_df <- human.med.utoronto.ca_complex_df[,c(2,1)]
write.table(human.med.utoronto.ca_complex_df, "/mnt/home3/reid/av638/tutorial/clustering/human.med.utoronto.ca_complex_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

# python benchmarkagain_complex.py dump.data.mci.I40 human.med.utoronto.ca_complex_df.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# python benchmarkagain_complex.py dump.data.mci.I40 human.med.utoronto.ca_complex_df.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'


# http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt
# cite:  https://www.embopress.org/doi/full/10.15252/msb.20167490

humap2_complexes_20200809 <- fread("/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809.txt")
humap2_complexes_20200809 <- data.frame(humap2_complexes_20200809)
humap2_complexes_20200809["genes"] <- unlist(lapply(strsplit(humap2_complexes_20200809$genenames, " "), function(x) paste0(x, collapse = ",")))
humap2_complexes_20200809_re <- humap2_complexes_20200809[,c(1,5)]
write.table(humap2_complexes_20200809_re, "/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809_re.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)


python3 benchmarkagain_complex.py dump.data.mci.I20 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
python3 benchmarkagain_complex.py dump.data.mci.I20 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'

# Relaxed cutoffs
python3 benchmarkagain_complex.py dump.data.mci.I20 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'

# ----> Note Taking column 4 from python3 benchmarkagain_complex.py   command output will help to understand KNOWN interactors 
# ----> Note Taking column 3 from python3 benchmarkagain_complex.py   command output will help to understand UNKNOWN interactors 

# ---------------------


summit_expdata_p1_exp <- read.table("/mnt/beegfs6/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/merged_out/intersected_one_all.txt", header = F, stringsAsFactors = F)
summit_expdata_p1_exp <- data.frame(summit_expdata_p1_exp[,c(1,2,5,6,7)])

colnames(summit_expdata_p1_exp) <- c("proteinA", "proteinB", "protein_Apeaks", "protein_Bpeaks", "Overlapspeaks")



summit_mydata_p1_exp <- rbind.data.frame(summit_expdata_p1_exp, summit_mydata_p1_cbh[,c(1:5)])
# -> 514+2 for test CAPs 
num_proteins = 516

summit_mydata_p1_exp["CAP_a"] <- summit_mydata_p1_exp$proteinA
summit_mydata_p1_exp["CAP_b"] <- summit_mydata_p1_exp$proteinB
summit_mydata_p1_exp_mg <- summit_mydata_p1_exp[,c(6,7,3:5)]
summit_mydata_p1_exp_mg["new_expected_value"] <- round((as.numeric(summit_mydata_p1_exp_mg$protein_Apeaks) * as.numeric(summit_mydata_p1_exp_mg$protein_Bpeaks))/250000) #lambda(expected_value)=(n*m)/N
summit_mydata_p1_exp_mg["pois_pval"] <- dpois(summit_mydata_p1_exp_mg$Overlapspeaks, summit_mydata_p1_exp_mg$new_expected_value)
#Perform Bonferroni Correction
summit_mydata_p1_exp_mg["Corrected_pois_pval"] <- p.adjust(summit_mydata_p1_exp_mg$pois_pval, method = "bonferroni", n = num_proteins * num_proteins)
summit_mydata_p1_exp_mg["Corrected_pois_pval"] <- ifelse(summit_mydata_p1_exp_mg$Corrected_pois_pval == 0, summit_mydata_p1_exp_mg$Corrected_pois_pval + 1e-300, summit_mydata_p1_exp_mg$Corrected_pois_pval)
summit_mydata_p1_exp_mg["Coassociation_score"] <- ifelse((summit_mydata_p1_exp_mg$Overlapspeaks > summit_mydata_p1_exp_mg$new_expected_value), -log(summit_mydata_p1_exp_mg$Corrected_pois_pval,10), log(summit_mydata_p1_exp_mg$Corrected_pois_pval,10))
summit_mydata_p1_exp_mg <- summit_mydata_p1_exp_mg %>% distinct()

summit_mydata_p1_exp_mg["proteinA"] <- summit_mydata_p1_exp_mg$CAP_a
summit_mydata_p1_exp_mg["proteinB"] <- summit_mydata_p1_exp_mg$CAP_b
summit_mydata_p1_exp_mg <- summit_mydata_p1_exp_mg[,c(10,11,3:9)]
#filter low count peaks
#sideA
summit_mydata_p1_exp_mg_filt <- summit_mydata_p1_exp_mg[which(summit_mydata_p1_exp_mg$protein_Apeaks > 5),]
#sideB
summit_mydata_p1_exp_mg_filt <- summit_mydata_p1_exp_mg_filt[which(summit_mydata_p1_exp_mg_filt$protein_Bpeaks > 5),]

#remove self-hits same protein on both sides
summit_mydata_p1_exp_mg_filt <- summit_mydata_p1_exp_mg_filt[which(summit_mydata_p1_exp_mg_filt$proteinA != summit_mydata_p1_exp_mg_filt$proteinB),]

#remove 0s coassociation score
summit_mydata_p1_exp_mg_filt2 <- summit_mydata_p1_exp_mg_filt[which(summit_mydata_p1_exp_mg_filt$Coassociation_score != 0),]

#remove +/- 1.4 coassociation
summit_mydata_p1_exp_mg_filt3 <- summit_mydata_p1_exp_mg_filt2[which(summit_mydata_p1_exp_mg_filt2$Coassociation_score > 1.4 | summit_mydata_p1_exp_mg_filt2$Coassociation_score < -1.4),]
#filter non-significant 0.05
# for (i in unique(summit_mydata_p1_exp_mg$proteinA)){
#   for (j in unique(summit_mydata_p1_exp_mg$proteinB)){
#     if (i != j){
#       print(paste0(i ,j))
#     }
#   }
# }

# write.table(summit_mydata_p1_exp_mg_filt3, "summit_mydata_p1_exp_mg_filt3.txt", sep="\t", append = F, quote = F, row.names = F)

summit_mydata_p1_exp_mg_re <- summit_mydata_p1_exp_mg_filt3[,c(1,2,9)]
colnames(summit_mydata_p1_exp_mg_re) <- c("protein1", "protein2", "score")
hist(summit_mydata_p1_exp_mg_re$score, breaks = 100)

write.table(summit_mydata_p1_exp_mg_re, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/summit_mydata_p1_exp_mg_re.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)


summit_mydata_p1_exp_mg_re.gt100 <- summit_mydata_p1_exp_mg_re[which(summit_mydata_p1_exp_mg_re$score >= 100),]
write.table(summit_mydata_p1_exp_mg_re.gt100, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/summit_mydata_p1_exp_mg_re.gt100.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)

summit_mydata_p1_exp_mg_re.gt10 <- summit_mydata_p1_exp_mg_re[which(summit_mydata_p1_exp_mg_re$score >= 10),]
write.table(summit_mydata_p1_exp_mg_re.gt10, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/summit_mydata_p1_exp_mg_re.gt10.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)


summit_mydata_p1_exp_mg_re_tab <- tidyr::pivot_wider(summit_mydata_p1_exp_mg_re, names_from = protein2, values_from = score)

#NA will come for those with no combinations matched
summit_mydata_p1_exp_mg_re_df <- as.data.frame(summit_mydata_p1_exp_mg_re_tab)
rownames(summit_mydata_p1_exp_mg_re_df) <- summit_mydata_p1_exp_mg_re_df$protein1
summit_mydata_p1_exp_mg_re_df <- summit_mydata_p1_exp_mg_re_df[,-1]
dim(summit_mydata_p1_exp_mg_re_df)

#replace NA with 0 to show the data was obsent
summit_mydata_p1_exp_mg_re_df[is.na(summit_mydata_p1_exp_mg_re_df)] <- 0

dim(summit_mydata_p1_exp_mg_re_df)

write.table(summit_mydata_p1_exp_mg_re_df, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/summit_mydata_p1_exp_mg_re_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

breaksListd <- seq(0, 1, by = 0.01)
library(pheatmap)
pheatmap_summit_mydata_p1_exp_mg_re_df <- pheatmap(summit_mydata_p1_exp_mg_re_df, 
                                                   na_col = "grey",
                                                   color = colorRampPalette(c("black", "white", "blue"))(length(breaksListd)),
                                                   clustering_distance_cols = "euclidean",
                                                   cluster_rows = T,
                                                   cluster_cols = T,
                                                   clustering_method = "ward.D")


# prepare hierarchical cluster
summit_mydata_p1_exp_mg_re_df_hc = hclust(dist(t(summit_mydata_p1_exp_mg_re_df)))
plot(summit_mydata_p1_exp_mg_re_df_hc)
# summit_mydata_p1_exp_mg_re_df_hc.pdf

# separate groups
summit_mydata_p1_exp_mg_re_df_cutree <- data.frame(groups=cutree(summit_mydata_p1_exp_mg_re_df_hc, k = 50))
summit_mydata_p1_exp_mg_re_df_cutree["caps"] <- rownames(summit_mydata_p1_exp_mg_re_df_cutree)
count_summit_mydata_p1_exp_mg_re_df_cutree <- count(summit_mydata_p1_exp_mg_re_df_cutree$groups)
colnames(count_summit_mydata_p1_exp_mg_re_df_cutree) <- c("complex", "num_of_caps")
plot(count_summit_mydata_p1_exp_mg_re_df_cutree$complex, count_summit_mydata_p1_exp_mg_re_df_cutree$num_of_caps)

summit_mydata_p1_exp_mg_re_df_cutree_merge <- ddply(summit_mydata_p1_exp_mg_re_df_cutree, .(groups), summarize, caps = toString(caps))
# summit_mydata_p1_exp_mg_re_df_cutree_merge[summit_mydata_p1_exp_mg_re_df_cutree_merge$caps %like% "CTCF", ]
summit_mydata_p1_exp_mg_re_df_cutree_merge[grepl("\\CTCF\\b", summit_mydata_p1_exp_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]
summit_mydata_p1_exp_mg_re_df_cutree_merge[grepl("\\SUZ12\\b", summit_mydata_p1_exp_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]
summit_mydata_p1_exp_mg_re_df_cutree_merge[grepl("\\NFYA\\b", summit_mydata_p1_exp_mg_re_df_cutree_merge$caps, ignore.case = TRUE),]

# use graph based clustering to identify complexes
# Cluster data with MCL
mcxload -abc summit_mydata_p1_exp_mg_re.gt100.txt --stream-mirror -write-tab data.tab -o data.mci 
# mcxload -abc summit_mydata_p1_exp_mg_re.gt10.txt --stream-mirror -write-tab data.tab -o data.mci 


mcl data.mci -I 1.4
mcl data.mci -I 2
mcl data.mci -I 4

mcxdump -icl out.data.mci.I14 -tabr data.tab -o dump.data.mci.I14
mcxdump -icl out.data.mci.I20 -tabr data.tab -o dump.data.mci.I20
mcxdump -icl out.data.mci.I40 -tabr data.tab -o dump.data.mci.I40

# Look at cluster sizes
cat dump.data.mci.I14  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I20  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I40  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'



# Obtain known protein complex data from  http://mips.helmholtz-muenchen.de/corum/#download

# Filter known complexes to only contain proteins from the ENCODE dataset
python3 complex_count_filter.py humanComplexes.txt summit_mydata_p1_exp_mg_re.txt > humanComplexes_filtered.txt

# Look at accuracy of predicted complexes
python3 mod_benchmark_complexes.py dump.data.mci.I14 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
python3 mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'

# Which complexes do we predict perfectly?
python3 mod_benchmark_complexes.py dump.data.mci.I20 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'
python3 mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'

# this script plots the senstivity and specificity to show improvements as inflation parameter increases
complex_analysis.R

library(readxl)
# http://human.med.utoronto.ca/php/data_download.php
# cite: https://www.cell.com/fulltext/S0092-8674(12)01006-9
human.med.utoronto.ca_TableS3 <- readxl::read_excel("/mnt/home3/reid/av638/tutorial/clustering/human.med.utoronto.ca_TableS3.xls", sheet = 1, col_names = FALSE, skip=1)
human.med.utoronto.ca_TableS3 <- data.frame(human.med.utoronto.ca_TableS3)
human.med.utoronto.ca_complex <- data.frame(complex=human.med.utoronto.ca_TableS3$...4)
human.med.utoronto.ca_complex <- data.frame(complex=human.med.utoronto.ca_complex[-1,])

human.med.utoronto.ca_complex_list <- lapply(strsplit(human.med.utoronto.ca_complex$complex, ","), function(x) paste0(gsub("_HUMAN","", x), collapse = ","))
human.med.utoronto.ca_complex_df <- do.call(rbind.data.frame, human.med.utoronto.ca_complex_list)
colnames(human.med.utoronto.ca_complex_df) <- "CAPs"
human.med.utoronto.ca_complex_df["complex_id"] <- paste0("complex_",rownames(human.med.utoronto.ca_complex_df))
human.med.utoronto.ca_complex_df <- human.med.utoronto.ca_complex_df[,c(2,1)]
write.table(human.med.utoronto.ca_complex_df, "/mnt/home3/reid/av638/tutorial/clustering/human.med.utoronto.ca_complex_df.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

# python benchmarkagain_complex.py dump.data.mci.I40 human.med.utoronto.ca_complex_df.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
# python benchmarkagain_complex.py dump.data.mci.I40 human.med.utoronto.ca_complex_df.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'


# http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt
# cite:  https://www.embopress.org/doi/full/10.15252/msb.20167490

humap2_complexes_20200809 <- fread("/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809.txt")
humap2_complexes_20200809 <- data.frame(humap2_complexes_20200809)
humap2_complexes_20200809["genes"] <- unlist(lapply(strsplit(humap2_complexes_20200809$genenames, " "), function(x) paste0(x, collapse = ",")))
humap2_complexes_20200809_re <- humap2_complexes_20200809[,c(1,5)]
write.table(humap2_complexes_20200809_re, "/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809_re.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)


python3 benchmarkagain_complex.py dump.data.mci.I20 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
python3 benchmarkagain_complex.py dump.data.mci.I20 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'

# Relaxed cutoffs
python3 benchmarkagain_complex.py dump.data.mci.I20 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.1'

# ----> Note Taking column 4 from python3 benchmarkagain_complex.py   command output will help to understand KNOWN interactors 
# ----> Note Taking column 3 from python3 benchmarkagain_complex.py   command output will help to understand UNKNOWN interactors 

