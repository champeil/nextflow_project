// This process is for fastqc for the fq.gz data
// author: laojp
// time: 2025.02.06
// position: SYSUCC bioinformatic platform
// usage
//      input: input_dir, input_suffix, output_suffix, sample_id, reads[list]

nextflow.enable.dsl=2

process process_fastp {
    label "fastp"
    tag "${sample_id}"
    publishDir path: "${params.result_dir}/fastp", mode: "copy"

    input:
        tuple val(input_dir), val(input_suffix), val(output_suffix), val(single_end), val(sample_id), path(reads)
    output:
        tuple val("${sample_id}"), 
            val("${single_end}"), 
            path("${sample_id}${output_suffix}*.fq.gz"), 
            path("${sample_id}${output_suffix}*json"), 
            path("${sample_id}${output_suffix}*html"), emit: fastp_result
    
    script:
        def fastp_cmd = ""
        if (single_end) {
            // single end
            fastp_cmd = """
                ${params.software_fastp}/fastp \
                    -i ${reads[0]} \
                    -o ${sample_id}${output_suffix}.fq.gz \
                    -j ${sample_id}${output_suffix}.json -h ${sample_id}${output_suffix}.html
            """
        } else {
            // paired end
            fastp_cmd = """
                ${params.software_fastp}/fastp \
                    -i ${reads[0]} -I ${reads[1]} \
                    -o ${sample_id}${output_suffix}_1.fq.gz -O ${sample_id}${output_suffix}_2.fq.gz \
                    -j ${sample_id}${output_suffix}.json -h ${sample_id}${output_suffix}.html
            """
        }

        // run fastp
        """
            #!/usr/bin/env bash
            ${fastp_cmd}
        """
}