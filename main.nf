#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import subworkflows
include { FASTA_PROCESSING  } from "$projectDir/subworkflows/fasta_processing/main"
include { XML_PROCESSING    } from "$projectDir/subworkflows/xml_processing/main"
include { RUN_BEAST         } from "$projectDir/subworkflows/run_beast/main"


// Disable resume when rerunning specific samples
if (params.rerun_beast_samples) {
    resume = false
}

workflow {

    // Set the workflow name
    def workflow_name = params.workflow_run_name ?: workflow.runName

    //
    log.info("""
    Running MitoDating. Workflow run name: $workflow_name
    """)

    // Define input channels
    Channel.fromPath(params.sample_metadata, checkIfExists: true).collect().set { metadata }
    Channel.fromPath(params.fasta, checkIfExists: true).set     { fasta }
    Channel.fromPath(params.priors, checkIfExists: true).collect().set    { priors }
    // Channel.fromPath(params.partition).set                      { partition }
    // Channel.fromPath(params.gff).set                            { gff }

    params.rerun_beast_samples = null  // e.g., "sample1,sample2,sample3" or null for full run

    if (params.rerun_beast_samples) {
        // Check if output directory exists
        if (!file(params.outdir).exists()) {
            error("Error: Output directory '${params.outdir}' not found. Run the full pipeline first.")
        }

        sample_list = params.rerun_beast_samples.split(',').collect { it.trim() }

        xml_files = Channel.fromPath("${params.outdir}/02_xml/*.xml")
            .filter { file ->
                sample_list.any { sample -> file.name.contains(sample) }
            }

        if (!xml_files) {
            error("No XML files found matching samples: ${params.rerun_beast_samples}")
        }

        RUN_BEAST(xml_files)
    } else {
        if ( params.fasta_processing.toBoolean() ) {
            FASTA_PROCESSING ( fasta, metadata )
        }

        if ( params.generate_XML.toBoolean() ) {
            XML_PROCESSING ( FASTA_PROCESSING.out.fastas, metadata, priors)
        }

        if ( params.run_beast.toBoolean() ) {
            RUN_BEAST( XML_PROCESSING.out.xml )
        }
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












