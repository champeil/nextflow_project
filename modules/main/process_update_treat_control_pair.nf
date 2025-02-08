// This process is for update the treat-control-pair sample with the raw treat-control-pair file
// author: laojp
// time: 2024.08.16
// position: SYSUCC bioinformatic platform
// usage
//      input: treat-normal-pair file in variant ${treat_normal_pair}
//      process: channel from each line of file with list [treat, control], 
//          then check the value channel with ",", if exist, then split and merge as the first element as final output respectively,
//          then collect the output files of each processed line which output to temp_file and combine as the new output
//      output: temp_file of each processed lines from raw treat-control-pair file

process update_treat_control_pair {
    label "update_treat_control_pair"
    tag "${treat_id}"

    input:
        tuple val(treat_id), val(control_id)

    output:
        path "single_sample_temp_file"
    
    script:
        """
            #!/usr/bin/env bash
            if [[ ${treat_id} == *","* ]]; then
                treat_id_merge=\$(echo \$(echo ${treat_id} | awk -v FS="," '{for(i=1;i<=NF;i++){print \$i}}' | head -n 1)_mer)
            else
                treat_id_merge=${treat_id}
            fi
            if [[ ${control_id} == *","* ]]; then
                control_id_merge=\$(echo \$(echo ${control_id} | awk -v FS="," '{for(i=1;i<=NF;i++){print \$i}}' | head -n 1)_mer)
            else
                control_id_merge=${control_id}
            fi
            echo -e "\${treat_id_merge}\t\${control_id_merge}" > "single_sample_temp_file"
        """

}
workflow workflow_update_treat_control_pair {
    raw_file=channel
        .fromPath("${params.result_dir}"+"/"+"${params.treat_control_pair}")
        .splitText()
        .map{ 
            it -> 
                def split_data = it.trim().split(/\s+/)
                return [split_data[0], split_data[1]]
        }
    update_treat_control_pair(raw_file) 
        .collectFile(name: params.treat_control_pair, storeDir: params.result_dir, newLine: false)
 }