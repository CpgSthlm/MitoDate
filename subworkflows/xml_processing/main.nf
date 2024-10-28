#! /usr/bin/env nextflow

include { GENERATEXML } from '../../modules/generatexml/main'

workflow XML_PROCESSING {
    take:
    fastas
    priors
    gff
    partition
    chainlength
    log_step
    partition_list
    nd_list
    taxon_set

    main:
    ch_versions = Channel.empty()

    GENERATEXML( fastas.flatten() , priors, gff, partition,
    chainlength, log_step, partition_list, nd_list, taxon_set
    )

    ch_versions = ch_versions.mix(GENERATEXML.out.versions)

    emit:
    xml             = GENERATEXML.out.xml
    versions        = GENERATEXML.out.versions
}

