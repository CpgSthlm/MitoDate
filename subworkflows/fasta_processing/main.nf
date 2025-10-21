#! /usr/bin/env nextflow

include { RENAMEFASTA } from '../../modules/rename_fastas/rename_fasta.nf'
include { SPLITFASTA } from '../../modules/splitfasta/splitfasta.nf'

workflow FASTA_PROCESSING {
    take:
    fasta
    metadata

    main:
    ch_versions                 = Channel.empty()

    // Run single sample dating
    if ( params.single_sample_dating ) {
        SPLITFASTA ( fasta, metadata )
        ch_versions             = ch_versions.mix(SPLITFASTA.out.versions)
    }
    else {
        RENAMEFASTA ( fasta, metadata )
        ch_versions                 = ch_versions.mix(RENAMEFASTA.out.versions)
    }


    emit:
    fastas      = params.single_sample_dating ? SPLITFASTA.out.fastas : RENAMEFASTA.out.renamed_fasta
    versions    = ch_versions
}



