#!/usr/bin/env nextflow

params.input_1 = "$projectDir/samplelist_1.txt"
params.outdir_1 = "$projectDir/allouts_1/"
params.input_2 = "$projectDir/samplelist_2.txt"
params.outdir_2 = "$projectDir/allouts_2/"
params.outdir = "$projectDir/merged_out/"
params.inputdir = "$projectDir/"
params.blacklist = "$projectDir/dummy_blacklist.bed"
params.finalout = "intersected_one_all.txt"
params.scripts = "$projectDir/"

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
    bedtools intersect -a ${bedfile_1}_filt_1.bed -b $blacklist -v >  ${bedfile_1}_filt_freeblack_1.bed
    """
}


process CONCATENATEBED_1 {
    cpus 8
    tag "Concat all files"
    publishDir params.outdir_1, mode: 'copy'

    input:
    path "*_filt_freeblack_1.bed"

    output:
    path "merged_bedfiles_1.bed"

    script:
    """
    cat  *_filt_freeblack_1.bed > merged_bedfiles_1.bed
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
    bedtools intersect -a ${bedfile_2}_filt_2.bed -b $blacklist -v >  ${bedfile_2}_filt_freeblack_2.bed
    """
}


process CONCATENATEBED_2 {
    cpus 8
    tag "Concat all files"
    publishDir params.outdir_2, mode: 'copy'

    input:
    path "*_filt_freeblack_2.bed"

    output:
    path "merged_bedfiles_2.bed"

    script:
    """
    cat *_filt_freeblack_2.bed > merged_bedfiles_2.bed
    """
}


process INTERSECTWITHMERGED{
    tag "Intersect one bed with other one files"
    publishDir params.outdir, mode: 'copy'
    input:
    path(pairs)
    path "merged_bedfiles_2.bed"

    output:
    path "${pairs[0]}_int_merged_bedfiles_2.bed"

    script:
    """
    bedtools intersect -a "${params.outdir_1}/${pairs[0]}"  -b "merged_bedfiles_2.bed" -wa -wb > "${pairs[0]}_int_merged_bedfiles_2.bed"
    """
}

process COUNTINTERSECTWITHMERGED{
    tag "Intersect one bed with other merged files"
    publishDir params.outdir, mode: 'copy'
    input:
    path(pairs)
    path "${pairs[0]}_int_merged_bedfiles_2.bed"

    output:
    stdout
    //nexflow will write standard output from here to separate file each time print in python executed

    script:
    """
    python ${params.scripts}paircounter.py --inputfile "${pairs[0]}_int_merged_bedfiles_2.bed" --currentpath ${params.inputdir} --processed_file_1 ${params.outdir_1} --processed_file_2 ${params.outdir_2} --processed_merged_2 ${params.outdir}
    """
}

process WRITEOUTPUT{
    tag "write count files to output"
    publishDir params.outdir, mode: 'copy'

    input:
    file(countfile)

    output:
    stdout

    script:
    """
    cat $countfile 
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
    cleanbedup_1_ch = FILTERBED_1(individual_beds_1_ch)
    blacklistfreebed_1_ch = REMOVEBLACKLISTBED_1(individual_beds_1_ch, cleanbedup_1_ch, params.blacklist)
    blacklistfreebed_1_ch
            .collect()
            .set { bedfile_1_ch } //for collecting all bedfiles name in one channel
    
    // bedfile_1_ch.view()

    concatbed_1_ch = CONCATENATEBED_1(bedfile_1_ch)
    // concatbed_1_ch
    //     .view()

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
    blacklistfreebed_2_ch = REMOVEBLACKLISTBED_2(individual_beds_2_ch, cleanbedup_2_ch, params.blacklist)
    blacklistfreebed_2_ch
            .collect()
            .set { bedfile_2_ch } //for collecting all bedfiles name in one channel
    
    // bedfile_2_ch.view()
    concatbed_2_ch = CONCATENATEBED_2(bedfile_2_ch)
    // concatbed_2_ch
    //     .view()

    // blacklistfreebed_1_ch.view()
    // blacklistfreebed_2_ch.view()

    pair_ch = blacklistfreebed_1_ch.combine(blacklistfreebed_2_ch)
    // pair_ch.view()

    pair_merged_ch = blacklistfreebed_1_ch.combine(blacklistfreebed_2_ch)
    // pair_merged_ch.view()

    intersectmerged_ch = INTERSECTWITHMERGED(pair_merged_ch, concatbed_2_ch)
    // intersectmerged_ch
    //     .view()
    countintersectmerged_ch = COUNTINTERSECTWITHMERGED(pair_merged_ch,intersectmerged_ch).collect()
    countintersectmerged_ch
        .collect()
        .set {outcountintersectmerged_ch}
    intersectout_ch = WRITEOUTPUT(outcountintersectmerged_ch)

}
