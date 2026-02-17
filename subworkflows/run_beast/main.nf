#! /usr/bin/env nextflow

include { BEAST             } from '../../modules/beast/main'
include { BEAST_LOG_PARSER  } from '../../modules/beastlogparse/main'
include { COLLECT_RESULTS   } from '../../modules/beastlogparse/collect_results'

workflow RUN_BEAST {
    take:
    xml

    main:
    ch_versions = Channel.empty()

    // Run BEAST
    BEAST( xml )
    ch_versions = ch_versions.mix(BEAST.out.versions)

    // Parse the BEAST logs
    BEAST_LOG_PARSER( BEAST.out.beast_out_log )
    ch_versions = ch_versions.mix(BEAST_LOG_PARSER.out.versions)




    // Compiling the results
    COLLECT_RESULTS( BEAST_LOG_PARSER.out.parsed_results.collect() )
    ch_versions = ch_versions.mix(COLLECT_RESULTS.out.versions)


    emit:
    beast_output        = BEAST.out.beast_out_log
    beast_trees         = BEAST.out.beast_out_trees
    parsed_results      = BEAST_LOG_PARSER.out.parsed_results
    parsed_log          = BEAST_LOG_PARSER.out.parsed_log
    age_summary         = COLLECT_RESULTS.out.age_summary
    versions            = ch_versions
}