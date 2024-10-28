#! /usr/bin/env nextflow

// Analysis script for MitoChronos pipeline
nextflow.enable.dsl = 2

// Import subworkflows

include { FASTA_PROCESSING  } from "$projectDir/subworkflows/split_fasta/main"
include { XML_PROCESSING    } from "$projectDir/subworkflows/xml_processing/main"
include { RUN_BEAST         } from "$projectDir/subworkflows/run_beast/main"



workflow {

    // Define workflow stages
    def recognized_workflow_stages = ['single_sample_dating','generate_XML','run_beast']

    // Check input
    def workflow_steps = params.steps.tokenize(",")
    if ( ! workflow_steps.every { it in recognized_workflow_stages } ) {
        error "Unrecognised workflow step in $params.steps ( $recognized_workflow_stages )"
    }

    //
    log.info("""
    Running MitoChronos.
    """)

    // Define input channels
    Channel.fromPath(params.fasta, checkIfExists: true).set     { fasta }
    Channel.fromPath(params.priors, checkIfExists: true).set    { priors }
    Channel.fromPath(params.partition).set                      { partition }
    Channel.fromPath(params.gff).set                            { gff }


    // Split fasta file for single sample dating

    if ( 'single_sample_dating' in workflow_steps ) {
        FASTA_PROCESSING ( fasta )
    }


    // Generate XML file
    if ( 'generate_XML' in workflow_steps ) {
        XML_PROCESSING ( FASTA_PROCESSING.out.fastas, priors, gff, partition,
        params.chainlength, params.log_step, params.partition_list, params.nd_list, params.taxon_set
        )
    }

    // Run BEAST
    if ( 'run_beast' in workflow_steps ) {
       RUN_BEAST( XML_PROCESSING.out.xml )
    }

}

// Workflow completion
workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Thank you for using MitoChronos.

        Results are located in the results folder.
        """)
    } else {
        log.info("""
        The pipeline completed unsuccessfully.

        Please read the error message. If you need help to solve your issue,
        feel free to reach out via slack or by opening an issue at
        .
        """)
    }
}












