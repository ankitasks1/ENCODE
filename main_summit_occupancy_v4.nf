#!/usr/bin/env nextflow

params.input_1 = "$projectDir/samplelist_1.txt"
params.outdir_1 = "$projectDir/allouts_1/"
params.outdir = "$projectDir/merged_out/"
params.inputdir = "$projectDir/"
params.blacklist_1 = "$projectDir/dummy_blacklist_1.bed"
params.gbinout = "summit_encode_out.txt"
params.scorefile = "total_sites_bound_summit_hsrfree_mat_st_h_re.txt"
params.scripts = "$projectDir/scripts/"
params.currentdir = "$projectDir/"
params.parameters = ""
params.help = ""
params.gbinsize = ""

// Show parameters message and exit
if (params.parameters){
    parameters()
    exit 0
}


if (params.help){
    helpmessage()
    exit 0
}

def parameters() {
log.info """\
    PAIRWISER ONE TO ALL - N F   P I P E L I N E
    =============================================
    Current Directory         : ${params.currentdir}
    Scripts Directory         : ${params.scripts}
    Samplesheet 1 Path        : ${params.input_1}
    Output 1 Directory        : ${params.outdir_1}
    Scripts location          : ${params.scripts}
    Output bin file           : ${params.gbinout}
    Output score file         : ${params.scorefile}
    """
    .stripIndent()
}

def helpmessage() {
log.info """\

usage: nextflow run main_summit_occupancy_v4.nf --input_1 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_occupancy/samplelist_1.txt --blacklist_1 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_occupancy/blacklistgrch38_ENCFF356LFX_1.bed --gbinsize 400
""".stripIndent()
}

process READSHEET_1{
    tag 'Reading content of sheet'
    publishDir params.outdir_1, mode: 'copy'
    input:
    path(beds)

    output:
    stdout

    script:
    """
    ls "$projectDir/$beds" | tr -d '\n'
    """
}

process FILTERBED_1{
    tag 'Filter peaks by -log10(q-value) > 1.3'
    publishDir params.outdir_1, mode: 'copy'
    input:
    path(bedfile_1)

    output:
    path "${bedfile_1}_filt_1.bed"

    script:
    """
    wc -l "${bedfile_1}"
    awk '{if(\$10 > 1.3) print \$0}'  ${bedfile_1} > ${bedfile_1}_filt_1.bed
    """
}

process REMOVEBLACKLISTBED_1{

    tag "Remove blacklisted regions from each bed files"
    publishDir params.outdir_1, mode: 'copy'
    input:
    path(bedfile_1)
    path "${bedfile_1}_filt_1.bed"
    path blacklist_1


    output:
    path "${bedfile_1}_filt_freeblack_1.bed"

    script:
    """
    bedtools intersect -a ${params.outdir_1}${bedfile_1}_filt_1.bed -b $blacklist_1 -v >  ${bedfile_1}_filt_freeblack_1.bed
    """
}

process GETSUMMIT_1{
    tag "get summit of peaks"
    publishDir params.outdir_1, mode: 'copy'
    
    input:
    path(bedfile_1)
    path "${bedfile_1}_filt_freeblack_1.bed"

    output:
    path "${bedfile_1}_summit_freeblack_sorted_merged_1.bed" //will send output files to work dir which can be detected by next processes. Should be only one output file is provided 

    script:
    """
    python3 ${params.scripts}get_summit_bed_occupancyspecific.py --inputfile "${bedfile_1}_filt_freeblack_1.bed" --processed_file ${params.outdir_1} --location processed_1
    sort -k1,1 -k2,2n "${params.outdir_1}${bedfile_1}_summit_freeblack_1.bed"  > "${params.outdir_1}${bedfile_1}_summit_freeblack_sorted_1.bed"
    bedtools merge -i "${params.outdir_1}${bedfile_1}_summit_freeblack_sorted_1.bed" > "${bedfile_1}_summit_freeblack_sorted_merged_1.bed"
    """
}

process COUNTFILTERED_1{
    tag "Count original and cleaned-up regions from each bed files"
    publishDir params.outdir_1, mode: 'copy'
    input:
    path(bedfile_1)
    path "${bedfile_1}_filt_freeblack_1.bed"

    output:
    stdout

    script:
    """
    python3 ${params.scripts}summit_peakcounter.py --inputfile "${bedfile_1}_summit_freeblack_sorted_merged_1.bed" --location other --processed_file ${params.outdir_1} --type cleaned_1
    python3 ${params.scripts}summit_peakcounter.py --inputfile "${bedfile_1}" --location current_1 --currentpath ${params.inputdir} --type original_1
    """
}

process MAKEQBINSFROMSUMMIT_1{
    tag "make genomic bins from each summit files"
    publishDir params.outdir_1, mode: 'copy'
    input:
    file "cleaned_1.count.txt"

    output:
    path "gbin_${params.gbinsize}_${params.gbinout}"

    script:
    """
    python3 ${params.scripts}makebins_from_summit_q_for_nf.py --inputfile "cleaned_1.count.txt" --gbinsize ${params.gbinsize} --output "gbin_${params.gbinsize}_${params.gbinout}" --path ${params.outdir_1}
    """
}



workflow{
    // Read samplesheet
    Channel
        .fromPath(params.input_1)
        .splitCsv(header: true)
        .map{ row -> file(row.beds) }
        .set{ sheet_1_ch }

    // Pass bed files one by one to channel
    individual_beds_1_ch = READSHEET_1(sheet_1_ch)

    // Filter bed files for quality
    cleanbedup_1_ch = FILTERBED_1(individual_beds_1_ch)


    // Remove blacklisted regions
    blacklistfreebed_1_ch = REMOVEBLACKLISTBED_1(individual_beds_1_ch, cleanbedup_1_ch, params.blacklist_1)
    blacklistfreebed_1_ch
            .collect()
            .set { bedfile_1_ch } //for collecting all bedfiles name in one channel

    // Create summit for each bed files
    get_summit_1_ch = GETSUMMIT_1(individual_beds_1_ch, bedfile_1_ch)
    get_summit_1_ch
            .collect()
            .set { summitbedfile_1_ch } //for collecting all bedfiles name in one channel
    

    // Count elements in filter bed files
    count_filtered_1_ch = COUNTFILTERED_1(individual_beds_1_ch, summitbedfile_1_ch)
    count_filtered_1_ch
        .collect()
        .set {count_collected_1_ch}

    //  Create specified bins
    makebins_ch = MAKEQBINSFROMSUMMIT_1(count_collected_1_ch)
    makebins_ch
        .collect()
        .set { makebinscollected_ch }

}
