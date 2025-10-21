#! /usr/bin/env nextflow

include { GENERATEXML } from '../../modules/generatexml/main'

workflow XML_PROCESSING {
    take:
    fasta
    metadata
    priors

    main:
    ch_versions = Channel.empty()

    GENERATEXML( fasta, metadata, priors )

    ch_versions = ch_versions.mix(GENERATEXML.out.versions)

    emit:
    xml             = GENERATEXML.out.xml
    versions        = GENERATEXML.out.versions
}

