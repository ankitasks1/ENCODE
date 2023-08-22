#!/usr/bin/env nextflow

params.input_1 = "$projectDir/samplelist_1.txt"
params.outdir_1 = "$projectDir/allouts_1/"
params.input_2 = "$projectDir/samplelist_2.txt"
params.outdir_2 = "$projectDir/allouts_2/"
params.outdir = "$projectDir/merged_out/"
params.inputdir = "$projectDir/"
params.blacklist = "$projectDir/dummy_blacklist.bed"
params.finalout = "intersected_one_all.txt"
params.intersectout = "intersected.counts.txt"
params.scripts = "$projectDir/"
params.currentdir = "$projectDir/"
params.parameters = ""
params.help = ""
params.hotspotregions = "$projectDir/merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50.txt"

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
    Samplesheet 2 Path        : ${params.input_2}
    Output 1 Directory        : ${params.outdir_1}
    Output 2 Directory        : ${params.outdir_2}
    Blacklist                 : ${params.blacklist}
    Blacklist                 : ${params.hotspotregions}
    """
    .stripIndent()
}

def helpmessage() {
log.info """\

usage: nextflow run main_summit_pairwise_one_to_all.nf --input_1 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/samplelist_1.txt --input_2 /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/samplelist_2.txt --blacklist /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/dummy_blacklist.bed --hotspotregions /Users/ankitverma/Documents/tutorial/nextflow/ENCODE/summit_pairwise/merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50.txt
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
    path blacklist


    output:
    path "${bedfile_1}_filt_freeblack_1.bed"

    script:
    """
    bedtools intersect -a ${params.outdir_1}${bedfile_1}_filt_1.bed -b $blacklist -v >  ${bedfile_1}_filt_freeblack_1.bed
    """
}

process REMOVEHOTSPOTREGIONS_1{

    tag "Remove hotspot regions from each bed files"
    publishDir params.outdir_1, mode: 'copy'
    input:
    path(bedfile_1)
    path "${bedfile_1}_filt_freeblack_1.bed"
    path hotspotregions


    output:
    path "${bedfile_1}_filt_freeblackhot_1.bed"

    script:
    """
    bedtools intersect -a ${params.outdir_1}${bedfile_1}_filt_freeblack_1.bed -b $hotspotregions -wa >  ${bedfile_1}_filt_freeblackhot_1.bed
    """
}

process GETSUMMIT_1{
    tag "get summit of peaks"
    publishDir params.outdir_1, mode: 'copy'
    
    input:
    path(bedfile_1)
    path "${bedfile_1}_filt_freeblackhot_1.bed"

    output:
    path "${bedfile_1}_summit_freeblackhot_sorted_merged_1.bed" //will send output files to work dir which can bed detected by nextf processes. Should be only one output file is provided 

    script:
    """
    python3 ${params.scripts}get_peaks_summit_bed.py --inputfile "${bedfile_1}_filt_freeblackhot_1.bed" --processed_file ${params.outdir_1} --location processed_1
    sort -k1,1 -k2,2n "${params.outdir_1}${bedfile_1}_summit_freeblackhot_1.bed"  > "${params.outdir_1}${bedfile_1}_summit_freeblackhot_sorted_1.bed"
    bedtools merge -i "${params.outdir_1}${bedfile_1}_summit_freeblackhot_sorted_1.bed" -c 12 -o distinct > "${bedfile_1}_summit_freeblackhot_sorted_merged_1.bed"
    """
}

process COUNTFILTERED_1{
    tag "Count original and cleaned-up regions from each bed files"
    publishDir params.outdir_1, mode: 'copy'
    input:
    path(bedfile_1)
    path "${bedfile_1}_filt_freeblackhot_1.bed"

    output:
    stdout

    script:
    """
    python3 ${params.scripts}peakcounter.py --inputfile "${bedfile_1}_summit_freeblackhot_sorted_merged_1.bed" --location other --processed_file ${params.outdir_1} --type cleaned_1
    python3 ${params.scripts}peakcounter.py --inputfile "${bedfile_1}" --location current_1 --currentpath ${params.inputdir} --type original_1
    """
}

process CONCATENATEBED_1 {
    cpus 8
    tag "Concat all files"
    publishDir params.outdir_1, mode: 'copy'

    input:
    path "*_filt_freeblackhot_1.bed"

    output:
    path "merged_bedfiles_1.bed"

    script:
    """
    cat  *_filt_freeblackhot_1.bed | sort -k1,1 -k2,2n > merged_bedfiles_1.bed
    """
}


process READSHEET_2{
    tag 'Reading content of sheet'
    publishDir params.outdir_2, mode: 'copy'
    input:
    path(beds)

    output:
    stdout

    script:
    """
    ls "$projectDir/$beds" | tr -d '\n'
    """
}

process FILTERBED_2{
    tag 'Filter peaks by -log10(q-value) > 1.3'
    publishDir params.outdir_2, mode: 'copy'
    input:
    path(bedfile_2)

    output:
    path "${bedfile_2}_filt_2.bed"

    script:
    """
    wc -l "${bedfile_2}"
    awk '{if(\$10 > 1.3) print \$0}' ${bedfile_2} > ${bedfile_2}_filt_2.bed
    """
}

process REMOVEBLACKLISTBED_2{
    tag "Remove blacklisted regions from each bed files"
    publishDir params.outdir_2, mode: 'copy'
    input:
    path(bedfile_2)
    path "${bedfile_2}_filt_2.bed"
    path blacklist


    output:
    path "${bedfile_2}_filt_freeblack_2.bed"

    script:
    """
    bedtools intersect -a ${params.outdir_2}${bedfile_2}_filt_2.bed -b $blacklist -v >  ${bedfile_2}_filt_freeblack_2.bed
    """
}

process REMOVEHOTSPOTREGIONS_2{

    tag "Remove hotspot regions from each bed files"
    publishDir params.outdir_2, mode: 'copy'
    input:
    path(bedfile_2)
    path "${bedfile_2}_filt_freeblack_2.bed"
    path hotspotregions


    output:
    path "${bedfile_2}_filt_freeblackhot_2.bed"

    script:
    """
    bedtools intersect -a ${params.outdir_2}${bedfile_2}_filt_freeblack_2.bed -b $hotspotregions -wa >  ${bedfile_2}_filt_freeblackhot_2.bed
    """
}

process GETSUMMIT_2{
    tag "get summit of peaks"
    publishDir params.outdir_2, mode: 'copy'
    
    input:
    path(bedfile_2)
    path "${bedfile_2}_filt_freeblackhot_2.bed"

    output:
    path "${bedfile_2}_summit_freeblackhot_sorted_merged_2.bed"

    script:
    """
    python3 ${params.scripts}get_peaks_summit_bed.py --inputfile "${bedfile_2}_filt_freeblackhot_2.bed" --processed_file ${params.outdir_2} --location processed_2
    sort -k1,1 -k2,2n "${params.outdir_2}${bedfile_2}_summit_freeblackhot_2.bed"  > "${params.outdir_2}${bedfile_2}_summit_freeblackhot_sorted_2.bed"
    bedtools merge -i "${params.outdir_2}${bedfile_2}_summit_freeblackhot_sorted_2.bed" -c 12 -o distinct > "${bedfile_2}_summit_freeblackhot_sorted_merged_2.bed"
    """
}

process COUNTFILTERED_2{
    tag "Count original and cleaned-up regions from each bed files"
    publishDir params.outdir_2, mode: 'copy'
    input:
    path(bedfile_2)
    path "${bedfile_2}_filt_freeblackhot_2.bed"

    output:
    stdout

    script:
    """
    python3 ${params.scripts}peakcounter.py --inputfile "${bedfile_2}_summit_freeblackhot_sorted_merged_2.bed" --location other --processed_file ${params.outdir_2} --type cleaned_2
    python3 ${params.scripts}peakcounter.py --inputfile "${bedfile_2}" --location current_2 --currentpath ${params.inputdir} --type original_2
    """
}

process CONCATENATEBED_2 {
    cpus 8
    tag "Concat all files"
    publishDir params.outdir_2, mode: 'copy'

    input:
    path "*_filt_freeblackhot_2.bed"

    output:
    path "merged_bedfiles_2.bed"

    script:
    """
    cat *_filt_freeblackhot_2.bed | sort -k1,1 -k2,2n > merged_bedfiles_2.bed
    """
}


process CLOSESTWITHMERGED{
    tag "Closest distance of one bed with other files"
    publishDir params.outdir, mode: 'copy'
    input:
    path(pairs)

    output:
    path "${pairs[0]}_int_merged_bedfiles_2.bed"

    script:
    """
    bedtools closest -a "${params.outdir_1}/${pairs[0]}"  -b "${params.outdir_2}/${pairs[1]}" -d > "${pairs[0]}_int_merged_bedfiles_2.bed"
    """
}


process COUNTINTERSECTWITHMERGED{
    tag "Intersect one bed with other merged files"
    publishDir params.outdir, mode: 'copy'
    input:
    path(pairs)
    path "${pairs[0]}_int_merged_bedfiles_2.bed"

    output:
    val "${params.intersectout}"

    script:
    """
    python3 ${params.scripts}distancecounter.py --inputfile "${pairs[0]}_int_merged_bedfiles_2.bed" --currentpath ${params.inputdir} --processed_file_1 ${params.outdir_1} --processed_file_2 ${params.outdir_2} --processed_merged_2 ${params.outdir} --outputbase ${params.intersectout}
    """
}

process WRITEOUTPUT{
    tag "write all count files to output"
    publishDir params.outdir, mode: 'copy'

    input:
    file 'combinecounts.py'

    output:
    val "${params.finalout}"

    script:
    """
    python3 ${params.scripts}combinecounts.py --currentpath ${params.inputdir} --processed_file_1 ${params.outdir_1} --processed_file_2 ${params.outdir_2} --processed_merged_2 ${params.outdir} --output ${params.finalout}
    """
}

workflow{

    Channel
        .fromPath(params.input_1)
        .splitCsv(header: true)
        .map{ row -> file(row.beds) }
        .set{ sheet_1_ch }

    // sheet_1_ch
    //     .view{ it }

    individual_beds_1_ch = READSHEET_1(sheet_1_ch)
    // individual_beds_1_ch
    //     .view()
    cleanbedup_1_ch = FILTERBED_1(individual_beds_1_ch)
    // cleanbedup_1_ch.view()
    blacklistfreebed_1_ch = REMOVEBLACKLISTBED_1(individual_beds_1_ch, cleanbedup_1_ch, params.blacklist)
    // blacklistfreebed_1_ch.view()

    hotlistfreebed_1_ch = REMOVEHOTSPOTREGIONS_1(individual_beds_1_ch, blacklistfreebed_1_ch, params.hotspotregions)
    // hotlistfreebed_1_ch.view()
    hotlistfreebed_1_ch
            .collect()
            .set { bedfile_1_ch } //for collecting all bedfiles name in one channel

    get_summit_1_ch = GETSUMMIT_1(individual_beds_1_ch, bedfile_1_ch)
    get_summit_1_ch
            .collect()
            .set { summitbedfile_1_ch } //for collecting all bedfiles name in one channel
    
    // summitbedfile_1_ch.view()
    count_filtered_1_ch = COUNTFILTERED_1(individual_beds_1_ch, summitbedfile_1_ch)
    // count_filtered_1_ch
    //     .view()

    concatbed_1_ch = CONCATENATEBED_1(summitbedfile_1_ch)
    concatbed_1_ch
        // .view()

    Channel
        .fromPath(params.input_2)
        .splitCsv(header: true)
        .map{ row -> file(row.beds) }
        .set{ sheet_2_ch }
    // sheet_2_ch
    //     .view{ it }

    individual_beds_2_ch = READSHEET_2(sheet_2_ch)
    // individual_beds_2_ch.view()
    cleanbedup_2_ch = FILTERBED_2(individual_beds_2_ch)
    // cleanbedup_2_ch.view()
    blacklistfreebed_2_ch = REMOVEBLACKLISTBED_2(individual_beds_2_ch, cleanbedup_2_ch, params.blacklist)
    // blacklistfreebed_2_ch.view()
    hotlistfreebed_2_ch = REMOVEHOTSPOTREGIONS_2(individual_beds_2_ch, blacklistfreebed_2_ch, params.hotspotregions)
    // hotlistfreebed_2_ch.view()
    hotlistfreebed_2_ch
            .collect()
            .set { bedfile_2_ch } //for collecting all bedfiles name in one channel

    get_summit_2_ch = GETSUMMIT_2(individual_beds_2_ch, bedfile_2_ch)
    get_summit_2_ch
            .collect()
            .set { summitbedfile_2_ch } //for collecting all bedfiles name in one channel

    count_filtered_2_ch = COUNTFILTERED_2(individual_beds_2_ch, summitbedfile_2_ch)
    // count_filtered_2_ch
    //     .view()

    // summitbedfile_2_ch.view()
    concatbed_2_ch = CONCATENATEBED_2(summitbedfile_2_ch)
    // concatbed_2_ch
    //     .view()

    // get_summit_1_ch.view()
    // get_summit_2_ch.view()

    pair_ch = get_summit_1_ch.combine(get_summit_2_ch)
    // pair_ch.view()

    pair_merged_ch = get_summit_1_ch.combine(concatbed_2_ch)
    pair_merged_ch.view() 

    closestmerged_ch = CLOSESTWITHMERGED(pair_merged_ch).collect()
    // closestmerged_ch
    //     .view()
   

    countintersectmerged_ch = COUNTINTERSECTWITHMERGED(pair_merged_ch,closestmerged_ch).collect()
    countintersectmerged_ch
        .collect()
        .set {outcountintersectmerged_ch}
    
    // // count_filtered_1_ch.view()
    // // count_filtered_2_ch.view()
    // // outcountintersectmerged_ch.view()

    outputfinal_ch = WRITEOUTPUT(outcountintersectmerged_ch)
}
