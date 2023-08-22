#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define input path
inputPath = "$baseDir/de.bam"

process indexBam {
    input:
    path bamPath

    script:
    """
    crun.samtools samtools index -b ${bamPath} -@ 16
    """
}

process convertToFastq {
    input:
    path bamPath

    output:
    path "${bamPath}.fastq"

    script:
    """
    crun.bamtools bamtools convert -in ${bamPath} -format fastq -out ${bamPath}.fastq
    """
}

process filterSeqtk {
    input:
    path fastqPath

    output:
    path "${fastqPath}_gte150bp.fastq"

    script:
    """
    crun.seqtk seqtk seq -L 150 ${fastqPath} > ${fastqPath}_gte150bp.fastq
    """
}

workflow {
    indexBam(inputPath)
    outfq=convertToFastq(inputPath)
    filterSeqtk(outfq)
}
