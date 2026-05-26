#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { RENAMEFASTA               } from "$projectDir/modules/rename_fastas/rename_fasta.nf"
include { SPLITFASTA                } from "$projectDir/modules/splitfasta/splitfasta.nf"
include { GENERATEXML               } from "$projectDir/modules/generatexml/main"
include { BEAST as BEAST_TIP_DATING } from "$projectDir/modules/beast/main"
include { BEAST_LOG_PARSER          } from "$projectDir/modules/beastlogparse/main"
include { COLLECT_RESULTS           } from "$projectDir/modules/beastlogparse/collect_results"
include { JOINT_TREE_META           } from "$projectDir/modules/joint_tree_meta/main"
include { JOINT_XML                 } from "$projectDir/modules/joint_xml/main"
include { BEAST as BEAST_JOINT_TREE } from "$projectDir/modules/beast/main"
include { TREE_ANNOTATOR            } from "$projectDir/modules/tree_annotator/main"
include { JOINT_TREE_VISUALIZE      } from "$projectDir/modules/joint_tree_visualize/main"


workflow {

    // Set the workflow name
    def workflow_name = params.workflow_run_name ?: workflow.runName

    // Validate run_mode parameter
    if (!(params.run_mode in ['tip_dating', 'rerun_samples', 'joint_tree'])) {
        error(
            "\tError: Invalid run_mode '${params.run_mode}'. Please set run_mode to 'tip_dating', 'rerun_samples', or 'joint_tree'."
        )
    }

    // Log the workflow run name and selected mode
    log.info("""\tWorkflow run name: $workflow_name""")
    log.info("""\tSelected run mode: $params.run_mode\n""")

    //////////////////////////////////////////////////////////////////////////////////////////////
    // MODE 1 — Full tip-dating workflow
    //////////////////////////////////////////////////////////////////////////////////////////////

    if (params.run_mode == 'tip_dating') {

        ch_metadata = Channel.fromPath(params.sample_metadata,  checkIfExists: true).collect()
        ch_fasta    = Channel.fromPath(params.fasta,            checkIfExists: true).collect()
        ch_priors   = Channel.fromPath(params.priors,           checkIfExists: true).collect()

        // FASTA PROCESSING
        if ( params.single_sample_dating.toBoolean() ) {
            SPLITFASTA( ch_fasta, ch_metadata )
            ch_fastas = SPLITFASTA.out.fastas.flatten()
        } else {
            RENAMEFASTA( ch_fasta, ch_metadata )
            ch_fastas = RENAMEFASTA.out.renamed_fasta
        }

        // XML GENERATION
        GENERATEXML( ch_fastas, ch_metadata, ch_priors )

        // RUN BEAST
        BEAST_TIP_DATING( GENERATEXML.out.xml )
        BEAST_LOG_PARSER( BEAST_TIP_DATING.out.beast_out_log )
        COLLECT_RESULTS( BEAST_LOG_PARSER.out.parsed_results.collect() )

    //////////////////////////////////////////////////////////////////////////////////////////////
    // MODE 2 — Rerun tip-dating for specific samples using existing FASTA files
    //////////////////////////////////////////////////////////////////////////////////////////////

    } else if (params.run_mode == 'rerun_samples') {

        log.info("""\tRerunning tip-dating for samples: ${params.rerun_samples}\n""")
        //log.info("""\tNote: Process cache is disabled for GENERATEXML, BEAST_TIP_DATING, BEAST_LOG_PARSER, and COLLECT_RESULTS\n""")

        if (!file(params.outdir).exists()) {
            error("\tError: Output directory '${params.outdir}' not found. Run the MODE 1 pipeline first.\n" +
                "\tMake sure to provide the same --outdir in the config file as in the original run," +
                "\tbecause the rerun will look for FASTA files in ${params.outdir}/01_fastas/\n")
        }

        ch_metadata = Channel.fromPath(params.sample_metadata, checkIfExists: true).collect()
        ch_priors   = Channel.fromPath(params.priors,          checkIfExists: true).collect()

        sample_list = params.rerun_samples.split(',').collect { it.trim() }

        fasta_files = Channel.fromPath("${params.outdir}/01_fastas/*.fasta")
            .filter { file ->
                sample_list.any { sample -> file.name.contains(sample) }
            }
            .ifEmpty { error("\tError: No FASTA files found in ${params.outdir}/01_fastas matching samples: ${params.rerun_samples}. Run MODE 1 with fasta processing first.") }

        // Regenerate XML, rerun BEAST, and collect results for the specified samples
        GENERATEXML( fasta_files, ch_metadata, ch_priors )
        BEAST_TIP_DATING( GENERATEXML.out.xml )
        BEAST_LOG_PARSER( BEAST_TIP_DATING.out.beast_out_log )
        COLLECT_RESULTS( BEAST_LOG_PARSER.out.parsed_results.collect() )

    //////////////////////////////////////////////////////////////////////////////////////////////
    // MODE 3 — Date calibrated Tree (after reviewing tip-dating results)
    //////////////////////////////////////////////////////////////////////////////////////////////

    } else if (params.run_mode == 'joint_tree') {

        log.info("""\tMaking date-calibrated tree using estimated tip dates from ${params.age_summary_file}\n""")

        if (!params.age_summary_file) {
            error(
                "\tError: Please provide --age_summary_file pointing to a Results_ageSummary.csv, e.g.:\n" +
                "\t--age_summary_file ${params.outdir}/05_age_summary/run_<timestamp>/Results_ageSummary.csv"
            )
        }
        if (!file(params.outdir).exists()) {
            error("\tError: Output directory '${params.outdir}' not found. Run the full pipeline first.")
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
        BEAST_JOINT_TREE( JOINT_XML.out.xml )

        // 4. Summarise the posterior tree sample into one MCC tree
        TREE_ANNOTATOR( BEAST_JOINT_TREE.out.beast_out_trees )

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

        Results are located in the folder: ${params.outdir}.
        """)
    } else {
        log.info("""
        The pipeline completed unsuccessfully.

        Please read the error message. If you need help to solve your issue,
        feel free to reach out by opening an issue at
        https://github.com/CpgSthlm/MitoDate/issues.
        """)
    }
}












