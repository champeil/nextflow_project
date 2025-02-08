// This process is for merge the replicate samples
// author: laojp
// time: 2024.08.15
// position: SYSUCC bioinformatic platform
// usage
//      input: raw treat-normal-pair in variant ${treat_normal_pair}
//      process: read file, then join each 
//      output: new merged fastq.gz file

process merge_raw_sample {
    label "merge_raw_sample"
    tag "${sample_id}"
    publishDir path: "${params.raw_fq_dir}", mode: "copy" 

    input:
        val sample_id

    output:
        path "*.fq.gz"
    when:
        sample_id =~ /,/

    script:
        """
            #!/usr/bin/env bash 
            sample_merge=\$(echo \$(echo ${sample_id} | awk -v FS="," '{for(i=1;i<=NF;i++){print \$i}}' | sort | head -n 1)_mer)
            if [[ ${params.single_end} == "true" ]]; then
                if [[ \$(find ${params.raw_fq_dir} -name \${sample_merge}*.fq.gz -type f) != 1 ]]; then
                    echo ${sample_id} | awk -v FS="," '{for(i=1;i<=NF;i++){print \$i}}' | while read id ; do 
                        zcat ${params.raw_fq_dir}/\${id}.fq.gz >> \${sample_merge}.fq.gz
                    done
                fi
            elif [[ ${params.single_end} == "false" ]]; then
                if [[ \$(find ${params.raw_fq_dir} -name \${sample_merge}*.fq.gz -type f) != 2 ]]; then
                    echo ${sample_id} | awk -v FS="," '{for(i=1;i<=NF;i++){print \$i}}' | while read id ; do 
                        cat ${params.raw_fq_dir}/\${id}_1.fq.gz >> \${sample_merge}_1.fq.gz
                        zcat ${params.raw_fq_dir}/\${id}_2.fq.gz >> \${sample_merge}_2.fq.gz
                    done
                fi
            fi
        """
}

workflow workflow_merge_raw_sample {
    raw_file=channel
        .fromPath("${params.result_dir}"+"/"+"${params.treat_control_pair}")
        .splitText()
        .map{ 
            it -> 
                def split_data = it.trim().split(/\s+/)
                return [split_data[0], split_data[1]]
        }
        .flatten()
    merge_raw_sample(raw_file)
}
