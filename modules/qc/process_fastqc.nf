#!/usr/bin/env nextflow
// This process is for fastqc for the fq.gz data
// author: laojp
// time: 2025.02.06
// position: SYSUCC bioinformatic platform
// usage
//      input: input_dir, input_suffix, output_suffix, sample_id, reads

nextflow.enable.dsl=2

process process_fastqc {
    label "fastqc"
    tag "${sample_id}"
    publishDir path: "${params.result_dir}/fastqc", mode: "copy"

    input:
        tuple val(input_dir), val(input_suffix), val(output_suffix), val(sample_id), path(reads)
    output:
        path "${sample_id}${output_suffix}*_fastqc.{zip,html}", emit: fastqc_report

    script:
        """
            #!/usr/bin/env bash
            for read in ${reads} ; do
                ${params.software_fastqc}/fastqc -t $task.cpus \${read} -o ./ 
            done
        """
}
