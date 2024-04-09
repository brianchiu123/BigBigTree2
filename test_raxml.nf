#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process TestRaxmlNG {

    tag "Raxml-NG Test"

    script:
    """
    sudo raxml-ng --help
    """
}

workflow{

    TestRaxmlNG()

}