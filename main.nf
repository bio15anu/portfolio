#!/usr/bin/env nextflow

inFasta = file(params.inFasta)

process oneline {

    publishDir "${params.wd}"

    input:
    file inFasta

    output:
    file "${params.outFasta}" into genomes

    """
    python3 $HOME/scripts/oneline.py ${params.inFasta} ${params.outFasta}
    """
}

process assembly_stats {

    input:
    file genome from genomes

    output:
    stdout result

    """
    python3 $HOME/scripts/assembly_stats.py ${genome}
    """
}

result.subscribe { println it.trim() }
