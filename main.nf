#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import subworkflows

include { FASTA_PROCESSING  } from "$projectDir/subworkflows/fasta_processing/main"
include { XML_PROCESSING    } from "$projectDir/subworkflows/xml_processing/main"
include { RUN_BEAST         } from "$projectDir/subworkflows/run_beast/main"



workflow {

    // Set the workflow name
    def workflow_name = params.workflow_run_name ?: workflow.runName

    //
    log.info("""
    Running MitoDating. Workflow run name: $workflow_name
    """)

    // Define input channels
    Channel.fromPath(params.sample_metadata, checkIfExists: true).set { metadata }
    Channel.fromPath(params.fasta, checkIfExists: true).set     { fasta }
    Channel.fromPath(params.priors, checkIfExists: true).set    { priors }
    // Channel.fromPath(params.partition).set                      { partition }
    // Channel.fromPath(params.gff).set                            { gff }

    if ( params.fasta_processing.toBoolean() ) {
        FASTA_PROCESSING ( fasta, metadata )
    }

    if ( params.generate_XML.toBoolean() ) {
        XML_PROCESSING ( FASTA_PROCESSING.out.fastas.collect(), metadata, priors)
    }

    if ( params.run_beast.toBoolean() ) {
        RUN_BEAST( XML_PROCESSING.out.xml )
    }

}

// Workflow completion
workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Thank you for using MitoDating.

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












