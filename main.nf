#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { RENAMEFASTA        } from "$projectDir/modules/rename_fastas/rename_fasta.nf"
include { SPLITFASTA         } from "$projectDir/modules/splitfasta/splitfasta.nf"
include { GENERATEXML        } from "$projectDir/modules/generatexml/main"
include { BEAST              } from "$projectDir/modules/beast/main"
include { BEAST_LOG_PARSER   } from "$projectDir/modules/beastlogparse/main"
include { COLLECT_RESULTS    } from "$projectDir/modules/beastlogparse/collect_results"


// Define parameters
params.rerun_beast_samples = null  // e.g., "sample1,sample2,sample3" or null for full run

// Disable resume when rerunning specific samples
if (params.rerun_beast_samples) {
    resume = false
}

workflow {

    // Set the workflow name
    def workflow_name = params.workflow_run_name ?: workflow.runName

    // Log the workflow name
    log.info("""
    Running MitoDate. Workflow run name: $workflow_name
    """)

    // Define input channels
    ch_metadata = Channel.fromPath(params.sample_metadata, checkIfExists: true).collect()
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
    ch_priors = Channel.fromPath(params.priors, checkIfExists: true).collect()

    //////////////////////////////////////////////////////////////////////////////////////////////

    if (params.rerun_beast_samples) {
        // Check if output directory exists
        if (!file(params.outdir).exists()) {
            error("Error: Output directory '${params.outdir}' not found. Run the full pipeline first.")
        }
        // Parse the sample list
        sample_list = params.rerun_beast_samples.split(',').collect { it.trim() }
        // Find XML files matching the sample list
        xml_files = Channel.fromPath("${params.outdir}/02_xml/*.xml")
            .filter { file ->
                sample_list.any { sample -> file.name.contains(sample) }
            }
        // Check if any XML files were found
        if (!xml_files) {
            error("No XML files found matching samples: ${params.rerun_beast_samples}")
        }
        // Run BEAST on selected samples
        BEAST( xml_files )
        BEAST_LOG_PARSER( BEAST.out.beast_out_log )
        COLLECT_RESULTS( BEAST_LOG_PARSER.out.parsed_results.collect() )
    } else {
        // FASTA PROCESSING
        if ( params.fasta_processing.toBoolean() ) {
            if ( params.single_sample_dating.toBoolean() ) {
                SPLITFASTA( ch_fasta, ch_metadata )
                ch_fastas = SPLITFASTA.out.fastas.flatten()
            } else {
                RENAMEFASTA( ch_fasta, ch_metadata )
                ch_fastas = RENAMEFASTA.out.renamed_fasta
            }
        }
        // XML GENERATION
        if ( params.generate_XML.toBoolean() ) {
            GENERATEXML( ch_fastas, ch_metadata, ch_priors )
        }

        // RUN BEAST
        if ( params.run_beast.toBoolean() ) {
            BEAST( GENERATEXML.out.xml )
            BEAST_LOG_PARSER( BEAST.out.beast_out_log )
            COLLECT_RESULTS( BEAST_LOG_PARSER.out.parsed_results.collect() )
        }
    }
}

// Workflow completion
workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Thank you for using MitoDate.

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












