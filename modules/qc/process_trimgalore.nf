// This process is for trimgalore to qc the data
// author: laojp
// time: 2024.08.20
// position: SYSUCC bioinformatic platform
// usage
//      input: treat-normal-pair file in variant ${treat_normal_pair}
//      process: channel from each line of file with list [treat, control], 
//          then split as the single file and then construct the trimgalore for each samples
//      output: trimgalore result

nextflow.enable.dsl=2

process trimgalore {
    label "trimgalore"
    tag "${sample_id}"
    publishDir path: "${params.result_dir}/trimgalore", mode: "copy"

    input:
        tuple val(sample_id), val(sample_suffix), val(input_dir)
    output:
        path "${sample_id}_trimgalore*fq.gz*"
    script:
        """
            #!/usr/bin/env bash
            if [ ! -d ${params.result_dir}/trimgalore ]; then
                mkdir ${params.result_dir}/trimgalore
            fi
            if [[ ${params.single_end} == "false" ]]; then
                if [[ \$(find ${params.result_dir}/trimgalore -name "${sample_id}_trimgalore*fq.gz" | wc -l) != 2 ]]; then
                    ${params.software_trimgalore} \
                        --paired \
                        -q 28 --phred33 --length 30 --stringency 3 --gzip --cores 10\
                        --path_to_cutadapt ${params.software_cutadapt} \
                        --basename ${sample_id}_trimgalore \
                        ${input_dir}/${sample_id}${sample_suffix}_1.fq.gz ${input_dir}/${sample_id}${sample_suffix}_2.fq.gz 
                    mv ${sample_id}${sample_suffix}_trimgalore_val_1.fq.gz ${sample_id}_trimgalore_1.fq.gz
                    mv ${sample_id}${sample_suffix}_trimgalore_val_2.fq.gz ${sample_id}_trimgalore_2.fq.gz
                fi
            elif [[ ${params.single_end} == "true" ]]; then 
                if [[ \$(find ${params.result_dir}/trimgalore -name "${sample_id}_trimgalore*fq.gz" | wc -l) != 1 ]]; then
                    ${params.software_trimgalore} \
                        -q 28 --phred33 --length 30 --stringency 3 --gzip --cores 10\
                        --path_to_cutadapt ${params.software_cutadapt} \
                        --basename ${sample_id}_trimgalore \
                        ${input_dir}/${sample_id}${sample_suffix}.fq.gz 
                    mv ${sample_id}${sample_suffix}_trimgalore_trimmed.fq.gz ${sample_id}_trimgalore.fq.gz
                fi
            fi
        """ 
}
workflow workflow_trimgalore {
    raw_file=channel
        .fromPath("${params.result_dir}/${params.treat_control_pair}")
        .splitText()
        .map{ 
            it -> 
                def split_data = it.trim().split(/\s+/)
                return [split_data[0], split_data[1]]
        }
        .flatten()
        .map{it -> return [it, "", "${params.result_dir}/${params.raw_fq_dir}"]}
    trimgalore(raw_file)
}