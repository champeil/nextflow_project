// This process is for markduplicate the bam file to remove duplication
// author: laojp
// time: 2024.08.22
// position: SYSUCC bioinformatic platform
// usage
//      input: treat-normal-pair file in variant ${treat_normal_pair}, and the input dir, and the suffix
//      process: channel from each line of file with list [treat, control], 
//          then split as the single file and then construct the markduplicate for each samples
//      output: markduplicate results

nextflow.enable.dsl=2

process markduplicate {
    label "markduplicate"
    tag "${sample_id}"
    publishDir path: "${params.result_dir}/markduplicate", mode: "copy"

    input:
        tuple val(sample_id), val(sample_suffix), val(input_dir)
    
    output:
        path "${sample_id}_marked.ba*"
    
    script:
        """
            #!/usr/bin/env bash
            if [ ! -d ${params.result_dir}/markduplicate ]; then
                mkdir -p ${params.result_dir}/markduplicate
            fi 
            if [[ \$(find ${params.result_dir}/markduplicate -name "${sample_id}_marked.bam" | wc -l) != 1 ]]; then
                ${params.software_gatk} MarkDuplicates \
                    -I ${input_dir}/${sample_id}${sample_suffix}.bam \
                    --REMOVE_DUPLICATES true \
                    --CREATE_INDEX true \
                    -O ${sample_id}_marked.bam \
                    -M ${sample_id}_marked.metrics
            fi
        """
}

workflow workflow_markduplicate {
    raw_file=channel
        .fromPath("${params.result_dir}/${params.treat_control_pair}")
        .splitText()
        .map{ 
            it -> 
                def split_data = it.trim().split(/\s+/)
                return [split_data[0], split_data[1]]
        }
        .flatten()
        .map{it -> return [it, "_bwa", "${params.result_dir}/bwa"]}
    markduplicate(raw_file)
}