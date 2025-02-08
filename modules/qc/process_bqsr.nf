// This process is for trimgalore to qc the data
// author: laojp
// time: 2024.08.20
// position: SYSUCC bioinformatic platform
// usage
//      input: treat-normal-pair file in variant ${treat_normal_pair}
//      process: channel from each line of file with list [treat, control], 
//          then split as the single file and then construct the bqsr for each samples
//      output: bqsr result

nextflow.enable.dsl=2

process bqsr {
    label "bqsr"
    tag "${sample_id}"
    publishDir path: "${params.result_dir}/bqsr", mode: "copy"

    input:
        tuple val(sample_id), val(sample_suffix), val(input_dir)
    output:
        path "${sample_id}_bqsr.ba*"
    script:
        """
            #!/usr/bin/env bash
            if [ ! -d ${params.result_dir}/bqsr ]; then
                mkdir -p ${params.result_dir}/bqsr
            fi
            if [[ \$(find ${params.result_dir}/bqsr -name "${sample_id}_bqsr.bam" | wc -l) != 1 ]]; then
                ${params.software_gatk} BaseRecalibrator \
                    -R ${params.reference_fasta} \
                    -I ${input_dir}/${sample_id}${sample_suffix}.bam \
                    --known-sites ${params.reference_dbsnp} \
                    --known-sites ${params.reference_gold_indel} \
                    -O ${sample_id}_recal.table
                ${params.software_gatk} ApplyBQSR \
                    -R ${params.reference_fasta} \
                    -I ${input_dir}/${sample_id}${sample_suffix}.bam \
                    -bqsr ${sample_id}_recal.table \
                    --create-output-bam-index true \
                    -O ${sample_id}_bqsr.bam
            fi
        """ 
}
workflow workflow_bqsr {
    raw_file=channel
        .fromPath("${params.result_dir}/${params.treat_control_pair}")
        .splitText()
        .map{ 
            it -> 
                def split_data = it.trim().split(/\s+/)
                return [split_data[0], split_data[1]]
        }
        .flatten()
        .map{it -> return [it, "_marked", "${params.result_dir}/markduplicate"]}
    bqsr(raw_file)
}