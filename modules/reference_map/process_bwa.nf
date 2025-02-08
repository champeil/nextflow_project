// This process is for mapping the qc-ed fq to reference with bwa
// author: laojp
// time: 2024.08.22
// position: SYSUCC bioinformatic platform
// usage
//      input: treat-normal-pair file in variant ${treat_normal_pair}, and the input dir, and the suffix
//      process: channel from each line of file with list [treat, control], 
//          then split as the single file and then construct the bwa for each samples
//      output: bwa results

nextflow.enable.dsl=2

process bwa {
    label "bwa"
    tag "${sample_id}"
    publishDir path: "${params.result_dir}/bwa", mode: "copy"

    input:
        tuple val(sample_id), val(sample_suffix), val(input_dir)
    
    output:
        path "${sample_id}_bwa.ba*"
    
    script:
        """
            #!/usr/bin/env bash
            if [ ! -d ${params.result_dir}/bwa ]; then
                mkdir -p ${params.result_dir}/bwa
            fi
            if [[ ${params.single_end} == "false" ]]; then
                if [[ \$(find ${params.result_dir}/bwa -name "\${sample_id}_bwa.bam" | wc -l) != 1 ]]; then
                    ${params.software_bwa} mem -M -t 16\
                        -R "@RG\\tID:${1}\\tSM:${1}\\tLB:WXS\\tPL:Illumina" \
                        ${params.reference_bwa} \
                        ${input_dir}/${sample_id}${sample_suffix}_1.fq.gz ${input_dir}/${sample_id}${sample_suffix}_2.fq.gz \
                        -o ${sample_id}_bwa.sam
                    ${params.software_samtools} sort ${sample_id}_bwa.sam -o ${sample_id}_bwa.bam -O bam
                    ${params.software_samtools} index ${sample_id}_bwa.bam
                fi
            elif [[ ${params.single_end} == "true" ]]; then 
                if [[ \$(find ${params.result_dir}/bwa -name "\${sample_id}_bwa.bam" | wc -l) != 1 ]]; then
                    ${params.software_bwa} mem -M -t 16\
                        -R "@RG\\tID:${1}\\tSM:${1}\\tLB:WXS\\tPL:Illumina" \
                        ${params.reference_bwa} \
                        ${input_dir}/${sample_id}${sample_suffix}.fq.gz \
                        -o ${sample_id}_bwa.sam
                    ${params.software_samtools} sort ${sample_id}_bwa.sam -o ${sample_id}_bwa.bam -O bam
                    ${params.software_samtools} index ${sample_id}_bwa.bam
                fi
            fi
        """
}


workflow workflow_bwa {
    raw_file=channel
        .fromPath("${params.result_dir}/${params.treat_control_pair}")
        .splitText()
        .map{ 
            it -> 
                def split_data = it.trim().split(/\s+/)
                return [split_data[0], split_data[1]]
        }
        .flatten()
        .map{it -> return [it, "_trimgalore", "${params.result_dir}/trimgalore"]}
    bwa(raw_file)
}