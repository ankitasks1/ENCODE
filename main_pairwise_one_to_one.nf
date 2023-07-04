#!/usr/bin/env nextflow

params.input_1 = "$projectDir/samplelist_1.txt"
params.outdir_1 = "$projectDir/allouts_1/"
params.input_2 = "$projectDir/samplelist_2.txt"
params.outdir_2 = "$projectDir/allouts_2/"
params.outdir = "$projectDir/merged_out/"
params.inputdir = "$projectDir/"
params.blacklist = "$projectDir/dummy_blacklist.bed"
params.finalout = "intersect_out.txt"

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
    path "${bedfile_1}_filt.bed"

    script:
    """
    wc -l "${bedfile_1}"
    awk '{if(\$10 > 1.3) print \$0}'  ${bedfile_1} > ${bedfile_1}_filt.bed
    """
}

process REMOVEBLACKLISTBED_1{

    tag "Remove blacklisted regions from each bed files"
    publishDir params.outdir_1, mode: 'copy'
    input:
    path(bedfile_1)
    path "${bedfile_1}_filt.bed"
    path blacklist


    output:
    path "${bedfile_1}_filt_freeblack.bed"

    script:
    """
    bedtools intersect -a ${bedfile_1}_filt.bed -b $blacklist -v >  ${bedfile_1}_filt_freeblack.bed
    """
}


process CONCATENATEBED_1 {
    cpus 8
    tag "Concat all files"
    publishDir params.outdir_1, mode: 'copy'
    input:
    path '*_filt_freeblack.bed'

    output:
    path 'merged_bedfiles.bed'

    script:
    """
    cat *_filt_freeblack.bed > merged_bedfiles.bed
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
    path "${bedfile_2}_filt.bed"

    script:
    """
    wc -l "${bedfile_2}"
    awk '{if(\$10 > 1.3) print \$0}' ${bedfile_2} > ${bedfile_2}_filt.bed
    """
}

process REMOVEBLACKLISTBED_2{
    tag "Remove blacklisted regions from each bed files"
    publishDir params.outdir_2, mode: 'copy'
    input:
    path(bedfile_2)
    path "${bedfile_2}_filt.bed"
    path blacklist


    output:
    path "${bedfile_2}_filt_freeblack.bed"

    script:
    """
    bedtools intersect -a ${bedfile_2}_filt.bed -b $blacklist -v >  ${bedfile_2}_filt_freeblack.bed
    """
}


process CONCATENATEBED_2 {
    cpus 8
    tag "Concat all files"
    publishDir params.outdir_2, mode: 'copy'
    input:
    path '*_filt_freeblack.bed'

    output:
    path 'merged_bedfiles.bed'

    script:
    """
    cat *_filt_freeblack.bed > merged_bedfiles.bed
    """
}


process INTERSECT{
    tag "Intersect one bed with other one files"
    publishDir params.outdir, mode: 'copy'
    input:
    path(pairs)

    output:
    stdout

    script:

    """
    BEDNAME_1=\$(basename "${pairs[0]}" "_peaks_id.bed_filt_freeblack.bed")
    BEDNAME_2=\$(basename "${pairs[1]}" "_peaks_id.bed_filt_freeblack.bed")
    BEDCLEANED_1="${pairs[0]}"
    BEDCLEANED_2="${pairs[1]}"
    COUNT_BEDCLEANED_1=\$(cat \${BEDCLEANED_1} | wc -l)
    COUNT_BEDCLEANED_2=\$(cat \${BEDCLEANED_2} | wc -l)
    BEDORG_1=\$(basename "${pairs[0]}" "_filt_freeblack.bed")
    BEDORG_2=\$(basename "${pairs[1]}" "_filt_freeblack.bed")
    COUNT_BEDORG_1=\$(cat "${params.inputdir}\${BEDORG_1}" | wc -l)
    COUNT_BEDORG_2=\$(cat "${params.inputdir}\${BEDORG_2}" | wc -l)
    INTERSECT=\$(bedtools intersect -a "${params.outdir_1}/${pairs[0]}" -b "${params.outdir_2}/${pairs[1]}" -wa -wb | cut -f4 | sort -k1,1 -u | wc -l)
    echo "\${BEDNAME_1}\t\${BEDNAME_2}\t\${COUNT_BEDORG_1}\t\${COUNT_BEDORG_2}\t\${COUNT_BEDCLEANED_1}\t\${COUNT_BEDCLEANED_2}\t\${INTERSECT}" 
    """
}


process WRITEOUTPUT{
    tag "Write output"
    publishDir params.outdir, mode: 'copy'

    input:
    file(intersect_out)

    output:
    path "${params.finalout}"

    script:
    """
    cat ${intersect_out} > ${params.finalout}
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
    // concatbed_1_ch.view()

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
    // concatbed_2_ch.view()

    // blacklistfreebed_1_ch.view()
    // blacklistfreebed_2_ch.view()

    pair_ch = blacklistfreebed_1_ch.combine(blacklistfreebed_2_ch)
    // pair_ch.view()

    intersect_ch = INTERSECT(pair_ch).collect()

    intersectout_ch = WRITEOUTPUT(intersect_ch)

}
