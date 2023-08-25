# ENCODE


# Pairwise clustering among CAPs
## Step1: Modify original bed file to add sample name inside bed file content
<code> python3 get_protein_info.py </code>

##### Note: this will generate new_script_to_add_id.sh that need to be executed

<code> ./new_script_to_add_id.sh </code>

## Select some chip-seq bed for unknown(test) purpose. Eg. CRAMP1
awk '{print $1"\t"$2"\t"$3"\t""peaks_"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t""CRAMP1_V5cutnrun"}' V5_CR1_S10_V5_S9_peaks.narrowPeak | grep chr | sort -k1,1 -k2,2n > CRAMP1_V5cutnrun_peaks_id.bed

## Create a bedfiles list
ls -1 *.bed > bedfiles.txt

## Run deduplication.py to concatenate multiple bedfile for same CAP
python3 deduplicate_bedfiles.py  bedfiles.txt samplelist_1.txt

## Since both the files are same, I just copy it
cp samplelist_1.txt samplelist_2.txt

## Step2: Run nextflow script (developed for nextflow version >23; java > 11) (Note: This nextflow script version 4 will generate upto summit file filtered by quality, blacklist and hotspot regions)

<code>  conda activate nextflow_v23  </code>

#### All files are provided in repostory except hotspot regions which can be obtained by transfer 

<code> /mnt/home3/reid/av638/ENCODE/slurm_sub_re.py nextflow run main_summit_pairwise_one_to_all_v4.nf --input_1 /mnt/home3/reid/av638/ENCODE/summit_pairwise/samplelist_1.txt --input_2 /mnt/home3/reid/av638/ENCODE/summit_pairwise/samplelist_2.txt --blacklist_1 /mnt/home3/reid/av638/ENCODE/summit_pairwise/blacklistgrch38_ENCFF356LFX_1.bed --blacklist_2 /mnt/home3/reid/av638/ENCODE/summit_pairwise/blacklistgrch38_ENCFF356LFX_2.bed --hotspotregions /mnt/home3/reid/av638/ENCODE/summit_pairwise/merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50.txt -c /mnt/home3/nextflow/gurdon.config  </code>


##### Note: It needs one python scripts: summit_peakcounter.py to be present in the same directory.

## Step 3: Perform closest distance shell script:
<code>   
./run_shell_script.sh
echo "Analysis Started"
./shell_summit_intersector.sh > shell_summit_qualblackhot1000_50_intersector_output.txt
echo "Analysis Finished"

##### can be submitted with slurm in a shell script eg. run_bin_peak.sh
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py ./run_shell_script.sh

</code>




## Step5: Run nextflow script (developed for nextflow version >23; java > 11) (Note: This nextflow script version 5 will generate upto end count file after filtered by quality, blacklist and hotspot regions)
<code> 
conda activate nextflow_v23
</code>

<code> 
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py nextflow run main_summit_pairwise_one_to_all_v5.nf --input_1 /mnt/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/samplelist_1.txt --input_2 /mnt/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/samplelist_2.txt --blacklist_1 /mnt/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/blacklistgrch38_ENCFF356LFX_1.bed --blacklist_2 /mnt/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/blacklistgrch38_ENCFF356LFX_2.bed --hotspotregions /mnt/home3/reid/av638/ENCODE/summit_pairwise/experiment_test/merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50.txt -c /mnt/home3/nextflow/gurdon.config 
</code>

##### Note: It needs two python scripts: summit_poolcounts.py, and summit_peakcounter.py to be present in the same directory.

## Step 6: Calculate coassociation score
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

summit_mydata_p1_cbh_mg_re <- summit_mydata_p1_cbh_mg_filt3[,c(1,2,9)]
colnames(summit_mydata_p1_cbh_mg_re) <- c("protein1", "protein2", "score")
hist(summit_mydata_p1_cbh_mg_re$score, breaks = 100)

write.table(summit_mydata_p1_cbh_mg_re, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/summit_mydata_p1_cbh_mg_re.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)


summit_mydata_p1_cbh_mg_re.gt100 <- summit_mydata_p1_cbh_mg_re[which(summit_mydata_p1_cbh_mg_re$score >= 100),]
write.table(summit_mydata_p1_cbh_mg_re.gt100, "/mnt/home3/reid/av638/ENCODE/summit_pairwise/summit_mydata_p1_cbh_mg_re.gt100.txt", sep="\t", append = F, quote = F, row.names = F, col.names = F)

#### A SUMMIT OVERLAP MATRIX can bed produced from here for heirarchical clustering


## Step7: Use graph based clustering to identify complexes
##### Cluster data with MCL
mcxload -abc summit_mydata_p1_cbh_mg_re.gt100.txt --stream-mirror -write-tab data.tab -o data.mci 

mcl data.mci -I 1.4
mcl data.mci -I 2
mcl data.mci -I 4

mcxdump -icl out.data.mci.I14 -tabr data.tab -o dump.data.mci.I14
mcxdump -icl out.data.mci.I20 -tabr data.tab -o dump.data.mci.I20
mcxdump -icl out.data.mci.I40 -tabr data.tab -o dump.data.mci.I40

#### Look at cluster sizes
cat dump.data.mci.I14  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I20  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'
cat dump.data.mci.I40  | perl -ne 'chomp;@a=split/\t/;print scalar @a, "\n"'

#### Obtain known protein complex data from  http://mips.helmholtz-muenchen.de/corum/#download

##### Filter known complexes to only contain proteins from the ENCODE dataset
python3 complex_count_filter.py humanComplexes.txt summit_mydata_p1_cbh_mg_re.txt > humanComplexes_filtered.txt

##### Look at accuracy of predicted complexes
python3 mod_benchmark_complexes.py dump.data.mci.I14 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
python3 mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'

##### Which complexes do we predict perfectly?
python3 mod_benchmark_complexes.py dump.data.mci.I20 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'
python3 mod_benchmark_complexes.py dump.data.mci.I40 humanComplexes_filtered.txt | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'

##### this script plots the senstivity and specificity to show improvements as inflation parameter increases
complex_analysis.R

library(readxl)

##### http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt
##### cite:  https://www.embopress.org/doi/full/10.15252/msb.20167490

humap2_complexes_20200809 <- fread("/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809.txt")
humap2_complexes_20200809 <- data.frame(humap2_complexes_20200809)
humap2_complexes_20200809["genes"] <- unlist(lapply(strsplit(humap2_complexes_20200809$genenames, " "), function(x) paste0(x, collapse = ",")))
humap2_complexes_20200809_re <- humap2_complexes_20200809[,c(1,5)]
write.table(humap2_complexes_20200809_re, "/mnt/home3/reid/av638/tutorial/clustering/humap2_complexes_20200809_re.txt", sep = "\t", quote = F, append = F, col.names = F, row.names = F)

python3 benchmarkagain_complex.py dump.data.mci.I20 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'
python3 benchmarkagain_complex.py dump.data.mci.I20 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 1 && $a[4] >= 1'

#### Relaxed cutoffs
python3 benchmarkagain_complex.py dump.data.mci.I20 humap2_complexes_20200809_re.txt  | perl -ne 'chomp;@a=split/\t/;print "$_\n" if $a[5] >= 0.5'

##### ----> Note Taking column 4 from python3 benchmarkagain_complex.py   command output will help to understand KNOWN interactors 
#### ----> Note Taking column 3 from python3 benchmarkagain_complex.py   command output will help to understand UNKNOWN interactors 



## Additional steps: Binding pattern of CAPs in genomic sites (To find highly occupied regions) OCCUPANCY MATRIX
### Hot spopt regions were also identified by using approach in which peak is divided by window size and rounded-off to defined as occupancy:
<code> python makebins_from_peaks_q.py --input fsample_list.txt --blacklist blacklistgrch38_ENCFF356LFX.bed --gbinsize 1000 --output gbin_1000_peaks_encode_out.txt </code>

### Hot spopt regions were also identified by using the following approach:

#### Run summit analysis with hotspot filtering
#### Search hotspots
#### Choose cutoff windows
#### Choose 1000 bp as suggested and select
#### ./find_hotspot_region.sh
bedtools makewindows -g hg38.chrom.sizes -w 1000 > hg38_1000bp.txt
awk '{print $1"\t"$2"\t"$3"\t"$1"%"$2"%"$3}' hg38_1000bp.txt >  hg38_1000bp_merged.txt
bedtools intersect -wa -wb -a MERGED_final_cat_peaking.sorted.bed -b hg38_1000bp_merged.txt > MERGED_final_cat_peaking.sorted_hg38_1000bp.bed
awk '{print $12"%"$16}' MERGED_final_cat_peaking.sorted_hg38_1000bp.bed  | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"%"$3"%"$4}' > MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound.bed

#### Remove histone and iva data
fgrep -f sample_to_removeinternal_from_list.txt MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound.bed -w -v > MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected.bed

#### Run this R script to count and merge the coordinates and filter less than 50 proteins inside (find_hotspot_bound_regions.R)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
setwd("/Users/ankitverma/Documents/Archivio2/ENCODE/K562/hotspotfree")
MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected <- read.table("/Users/ankitverma/Documents/Archivio2/ENCODE/K562/hotspotfree/MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected.bed", header = F, stringsAsFactors = F)
detach("package:dplyr")
library(plyr)
head(MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected)
count_merged_final_cat_peaking.sorted_hg38_1000bp_justbound_selected <- count(MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected, "V2")

head(MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected)
clustered_MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected <- ddply(MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected, .(V2), summarize, Proteins = toString(V1))
#Identify same regions bound by multiple proteins
clustered_MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected$Proteins <- gsub(" ", "", clustered_MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected$Proteins)
dim(clustered_MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected)
head(clustered_MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected)

merged_encode_peak_count_1000bp_clustered_count <- merge(clustered_MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected, count_merged_final_cat_peaking.sorted_hg38_1000bp_justbound_selected, by="V2")
colnames(merged_encode_peak_count_1000bp_clustered_count) <- c("row", "Proteins", "freq")
library(splitstackshape)
merged_encode_peak_count_1000bp_clustered_count_sort <- cSplit(merged_encode_peak_count_1000bp_clustered_count, "row", "%")
head(merged_encode_peak_count_1000bp_clustered_count_sort)
merged_encode_peak_count_1000bp_clustered_count_sort <- merged_encode_peak_count_1000bp_clustered_count_sort[,c(3:5,1,2)]
merged_encode_peak_count_1000bp_clustered_count_sort <- merged_encode_peak_count_1000bp_clustered_count_sort[order(merged_encode_peak_count_1000bp_clustered_count_sort$row_1, merged_encode_peak_count_1000bp_clustered_count_sort$row_2),]
merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50 <- merged_encode_peak_count_1000bp_clustered_count_sort[which(merged_encode_peak_count_1000bp_clustered_count_sort$freq <= 50),]
head(merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50)
dim(merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50)
write.table(merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50, "merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)


Note: merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50.txt contains regions which are 1000 bp and bound by >=50 CAPs. So this file will be used for selecting only those peaks from all CAPs which intersect with these regions.




