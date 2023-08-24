#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.inputPath = null
params.onOmics = null
params.samToolsContainer = null
params.samToolsCpu = null
params.samToolsMemory = null
params.bamToolsContainer = null
params.bamToolsCpu = null
params.bamToolsMemory = null
params.seqtkContainer = null
params.seqtkCpu = null
params.seqtkMemory = null

publish_dir = { params.onOmics ? "/mnt/workflow/pubdir" : "." }

process indexBam {
    container params.samToolsContainer
    cpus params.samToolsCpu
    memory params.samToolsMemory
    publishDir publish_dir
    
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
    container params.bamToolsContainer
    cpus params.bamToolsCpu
    memory params.bamToolsMemory
    publishDir publish_dir
    
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
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir
    
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
