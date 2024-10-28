#! /usr/bin/env nextflow

include { BEAST } from '../../modules/beast/main'

workflow RUN_BEAST {
    take:
    xml

    main:
    ch_versions = Channel.empty()

    BEAST( xml.flatten() )

    ch_versions = ch_versions.mix(BEAST.out.versions)

    emit:
    beast_output        = BEAST.out.beast_output
    versions            = BEAST.out.versions
}