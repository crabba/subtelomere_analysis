#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
params.inputPath

process indexBam {
    container 'staphb/samtools:latest'

    publishDir "/mnt/workflow/pubdir"

    input:
        path bamPath

    output:
	    path bamPath

    script:
    """
    samtools index -b ${bamPath} -@ 16
    """
}

process convertToFastq {
    container 'biocontainers/bamtools:2.4.0'

    publishDir "/mnt/workflow/pubdir"

    input:
        path bamPath

    output:
        path "${bamPath}.fastq"

    script:
    """
    bamtools convert -in ${bamPath} -format fastq -out ${bamPath}.fastq
    """
}

process filterSeqtk {
    container 'staphb/seqtk:latest'

    publishDir "/mnt/workflow/pubdir"
    
    input:
        path fastqPath

    output:
        path "${fastqPath}_gte150bp.fastq"

    script:
    """
    seqtk seq -L 150 ${fastqPath} > ${fastqPath}_gte150bp.fastq
    """
}

workflow {
    def input_ch = Channel.fromPath(params.inputPath)
    indexBam(input_ch) | convertToFastq | filterSeqtk
}
