#! /usr/bin/env nextflow

include { SPLITFASTA } from '../../modules/splitfasta/main'

workflow FASTA_PROCESSING {
    take:
    fasta

    main:
    ch_versions                 = Channel.empty()

    SPLITFASTA ( fasta )
    ch_versions                 = ch_versions.mix(SPLITFASTA.out.versions)

    emit:
    fastas                      = SPLITFASTA.out.fastas
}



