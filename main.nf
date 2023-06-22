#!/usr/bin/env nextflow

params.input = "$projectDir/fsample_list.txt"
params.output = ""
params.outdir = "$projectDir/allouts"
params.blacklist = "$projectDir/dummy_blacklist.bed"
params.gbinsize = 1000
params.evalue = 1.3
params.peakscore = 0
params.scorebased = ""
params.path = "$projectDir"
params.pyscript = "$projectDir/makebins_from_peaks_q.py"



log.info """\
    P E A K B I N M A R K E R - P I P E L I N E
    =========================================
    ChIp-Seq file Path        : ${params.input}
    Blacklisted regions       : ${params.blacklist}    
    Output file name          : ${params.output}
    Output Files Directory    : ${params.outdir}
    Pyhton script to execute  : ${params.pyscript}
    Genomic Bin size.         : ${params.gbinsize}
    E-value-based analy       : ${params.evalue}
    """
    .stripIndent()

process PEAKBINMARKER {
    publishDir params.outdir, mode: 'copy'
    debug true
    input:
    path input
    path blacklist
    path pyscript
    val gbinsize

    output:
    path "${params.output}"

    script:
    """
    rm -r ${params.outdir}
    mkdir ${params.outdir}
    echo "$PWD"
    echo "${params.input}" > ${params.output}
    python $pyscript --input $input --filepath /workspace/gitpod/nf-training/data/peakbinmarker/ --blacklist $blacklist --gbinsize $gbinsize --output ${params.output}
    """
}

workflow {
    output_ch = PEAKBINMARKER(params.input, params.blacklist, params.pyscript, params.gbinsize)
    output_ch
}


