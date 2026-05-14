#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { RENAMEFASTA               } from "$projectDir/modules/rename_fastas/rename_fasta.nf"
include { SPLITFASTA                } from "$projectDir/modules/splitfasta/splitfasta.nf"
include { GENERATEXML               } from "$projectDir/modules/generatexml/main"
include { BEAST as BEAST_DATING     } from "$projectDir/modules/beast/main"
include { BEAST_LOG_PARSER          } from "$projectDir/modules/beastlogparse/main"
include { COLLECT_RESULTS           } from "$projectDir/modules/beastlogparse/collect_results"
include { JOINT_TREE_META           } from "$projectDir/modules/joint_tree_meta/main"
include { JOINT_XML                 } from "$projectDir/modules/joint_xml/main"
include { BEAST as BEAST_TREE       } from "$projectDir/modules/beast/main"
include { TREE_ANNOTATOR            } from "$projectDir/modules/tree_annotator/main"
include { JOINT_TREE_VISUALIZE      } from "$projectDir/modules/joint_tree_visualize/main"


// Define parameters
params.rerun_beast_samples = null  // e.g., "sample1,sample2,sample3" or null for full run
params.final_dated_tree    = false // run the joint dated-tree step as a standalone mode
params.age_summary_file    = null  // required with --final_dated_tree: path to Results_ageSummary.csv

// Disable resume when rerunning specific samples
if (params.rerun_beast_samples) {
    resume = false
}

workflow {

    // Set the workflow name
    def workflow_name = params.workflow_run_name ?: workflow.runName

    log.info("""
    Running MitoDate. Workflow run name: $workflow_name
    """)

    //////////////////////////////////////////////////////////////////////////////////////////////
    // MODE 1 — Main run of the pipeline
    //////////////////////////////////////////////////////////////////////////////////////////////

    if (!params.rerun_beast_samples && !params.final_dated_tree.toBoolean()) {
        ch_metadata = Channel.fromPath(params.sample_metadata,  checkIfExists: true).collect()
        ch_fasta    = Channel.fromPath(params.fasta,            checkIfExists: true)
        ch_priors   = Channel.fromPath(params.priors,           checkIfExists: true).collect()

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
            BEAST_DATING( GENERATEXML.out.xml )
            BEAST_LOG_PARSER( BEAST_DATING.out.beast_out_log )
            COLLECT_RESULTS( BEAST_LOG_PARSER.out.parsed_results.collect() )
        }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // MODE 2 — Rerun BEAST for specific samples
    //////////////////////////////////////////////////////////////////////////////////////////////

    } else if (params.rerun_beast_samples) {
        if (!file(params.outdir).exists()) {
            error("Output directory '${params.outdir}' not found. Run the MODE 1 pipeline first.")
        }
        sample_list = params.rerun_beast_samples.split(',').collect { it.trim() }
        xml_files = Channel.fromPath("${params.outdir}/02_xml/*.xml")
            .filter { file ->
                sample_list.any { sample -> file.name.contains(sample) }
            }
            .ifEmpty { error("No XML files found matching samples: ${params.rerun_beast_samples}") }
        BEAST_DATING( xml_files )
        BEAST_LOG_PARSER( BEAST_DATING.out.beast_out_log )
        COLLECT_RESULTS( BEAST_LOG_PARSER.out.parsed_results.collect() )

    //////////////////////////////////////////////////////////////////////////////////////////////
    // MODE 3 — Date calibrated Tree (after reviewing tip-dating results)
    //////////////////////////////////////////////////////////////////////////////////////////////

    } else if (params.final_dated_tree.toBoolean()) {
        if (!params.age_summary_file) {
            error(
                "Please provide --age_summary_file pointing to a Results_ageSummary.csv from a completed run, e.g.:\n" +
                "--age_summary_file ${params.outdir}/05_age_summary/run_<timestamp>/Results_ageSummary.csv"
            )
        }
        if (!file(params.outdir).exists()) {
            error("Output directory '${params.outdir}' not found. Run the full pipeline first.")
        }

        ch_metadata    = Channel.fromPath(params.sample_metadata, checkIfExists: true).collect()
        ch_age_summary = Channel.fromPath(params.age_summary_file, checkIfExists: true)
        ch_fasta       = Channel.fromPath(params.fasta,            checkIfExists: true)

        // 1. Merge estimated dates into metadata and rename FASTA accordingly
        JOINT_TREE_META( ch_age_summary, ch_metadata, ch_fasta )

        // 2. Generate BEAST XML — no priors table needed (all samples are now dated)
        JOINT_XML(
            JOINT_TREE_META.out.joint_fasta,
            JOINT_TREE_META.out.joint_meta
        )

        // 3. Run BEAST for the joint dated tree
        BEAST_TREE( JOINT_XML.out.xml )

        // 4. Summarise the posterior tree sample into one MCC tree
        TREE_ANNOTATOR( BEAST_TREE.out.beast_out_trees )

        // 5. Plot the annotated tree coloured by Group-By
        JOINT_TREE_VISUALIZE(
            TREE_ANNOTATOR.out.annotated_tree,
            JOINT_TREE_META.out.joint_meta
        )
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
        https://github.com/CpgSthlm/MitoDate/issues.
        """)
    }
}












