#! /usr/bin/env nextflow

include { BEAST             } from '../../modules/beast/main'
include {BEAST_LOG_PARSER   } from '../../modules/beastlogparse/main'

workflow RUN_BEAST {
    take:
    xml

    main:
    ch_versions = Channel.empty()

    BEAST( xml )
    ch_versions = ch_versions.mix(BEAST.out.versions)

    BEAST_LOG_PARSER( BEAST.out.beast_out_log )
    ch_versions = ch_versions.mix(BEAST_LOG_PARSER.out.versions)


    emit:
    beast_output        = BEAST.out.beast_out_log
    beast_trees         = BEAST.out.beast_out_trees
    versions            = ch_versions
}