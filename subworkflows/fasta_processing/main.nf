#! /usr/bin/env nextflow

include { RENAMEFASTA } from '../../modules/rename_fastas/rename_fasta.nf'

workflow FASTA_PROCESSING {
    take:
    fasta
    metadata

    main:
    ch_versions                 = Channel.empty()

    RENAMEFASTA ( fasta, metadata )
    ch_versions                 = ch_versions.mix(RENAMEFASTA.out.versions)


    emit:
    renamed_fasta                      = RENAMEFASTA.out.renamed_fasta
    versions                           = ch_versions
}



