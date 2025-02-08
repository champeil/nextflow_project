#!/usr/bin/env nextflow
// this script is for running rnaseq pipeline
// author: laojp
// time: 2025.02.05
// position: SYSUCC bioinformatic platform
// usage: 
//  nextflow run rnaseq.nf -c RNAseq.config 
//   -with-report report.html -with-trace -with-timeline timeline.html -with-dag dag.html

// base dir configure
nextflow.enable.dsl=2

// Use baseDir for relative paths
include { process_fastqc; process_fastqc as process_fastqc_clean } from "${projectDir}/../../modules/qc/process_fastqc.nf"
include { process_multiqc; process_multiqc as process_multiqc_clean } from "${projectDir}/../../modules/qc/process_multiqc.nf"
include { process_fastp } from "${projectDir}/../../modules/qc/process_fastp.nf"
include { process_star } from "${projectDir}/../../modules/reference_map/process_star.nf"
include { process_featurecount; process_featurecount_merge } from "${projectDir}/../../modules/quantify/process_featurecount.nf"

// read channel
Channel
    .fromPath("${params.sample_list}")
    .splitCsv(sep: ',', header: true)
    .map { row -> 
        if(row.read2 == 'false'){
            return tuple(row.sample_id, [ file("${row.read1}") ], true)
        } else {
            return tuple(row.sample_id, [ file("${row.read1}"), file("${row.read2}") ], false)
        }}
    .set { sample_list }

// workflow
workflow {
    // fastqc
    fastqc_result = sample_list
        .map { read_input -> tuple("", "", "", read_input[0], read_input[1]) }
        | process_fastqc

    // multiqc
    fastqc_result
        .collect()
        .map { fastqc_output -> tuple(fastqc_output, "", "_raw") }
        | process_multiqc

    // fastp
    fastp_result = sample_list
        .map { read_input -> tuple("", "", "_fastp", read_input[2], read_input[0], read_input[1]) }
        | process_fastp

    // fastqc
    fastqc_result_fastp = fastp_result
        .map { sample_id, single_end, clean_reads, json, html ->
            tuple("${params.result_dir}/fastp", "", "", sample_id, clean_reads) }
        | process_fastqc_clean

    // multiqc
    fastqc_result_fastp
        .collect()
        .map { fastqc_output -> tuple(fastqc_output, "", "_fastp") }
        | process_multiqc_clean

    // STAR
    star_result = fastp_result
        .map { sample_id, single_end, clean_reads, json, html ->
            tuple("${params.result_dir}/fastp", "_fastp", "_star", single_end, sample_id, clean_reads) }
        | process_star

    // featurecount
    featurecount_result = star_result
        .map { sample_id, single_end, gene_bam, transcript_bam, map_log, log ->
            tuple("${params.result_dir}/star", "_star", "_featurecount", single_end, sample_id, gene_bam) }
        | process_featurecount

    featurecount_result.featurecount_output
        .collect()
        .map { featurecount_res -> 
            tuple("${params.result_dir}/featurecount", "_featurecount", "_merge", featurecount_res)
        }
        | process_featurecount_merge
    
}
