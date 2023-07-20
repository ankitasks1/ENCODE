# ENCODE


# Pairwise clustering among CAPs
## Step1: Modify original bed file to add sample name inside bed file content
<code> python3 get_protein_info.py </code>
<code> ./new_script_to_add_id.sh </code>

## Step2: Run nextflow script (developed for nextflow version >23; java > 11)
<code> /mnt/home3/reid/av638/ENCODE/slurm_sub_re.py nextflow run main_pairwise_one_to_all.nf --input_1 /mnt/home3/reid/av638/ENCODE/pairwise/samplelist_1.txt --input_2 /mnt/home3/reid/av638/ENCODE/pairwise/samplelist_2.txt --blacklist /mnt/home3/reid/av638/ENCODE/pairwise/blacklistgrch38_ENCFF356LFX.bed -c /mnt/home3/nextflow/gurdon.config  </code>


#### Note: It needs three python scripts: paircounter.py, combinecounts.py and peakcounter.py to be present in the same directory.

# Binding pattern of CAPs in genomic sites
<code> python makebins_from_peaks_q.py --input fsample_list.txt --blacklist blacklistgrch38_ENCFF356LFX.bed --gbinsize 1000 --output gbin_1000_peaks_encode_out.txt </code>

##### can be submitted with slurm in a shell script eg. run_bin_peak.sh
