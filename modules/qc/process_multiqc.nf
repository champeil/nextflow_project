#!/usr/bin/env nextflow
// This process is for multiqc for the fastqc data
// author: laojp
// time: 2025.02.06
// position: SYSUCC bioinformatic platform
// usage
//      input: fastqc_report, input_suffix, output_suffix

nextflow.enable.dsl=2

process process_multiqc {
    label "multiqc"
    tag "multiqc"
    publishDir path: "${params.result_dir}/multiqc", mode: "copy"

    input:
        tuple path(fastqc_report), val(input_suffix), val(output_suffix)
    output:
        path "multiqc${output_suffix}*"
    
    script:
        """
            #!/usr/bin/env bash
            ${params.software_multiqc}/multiqc -n "multiqc${output_suffix}" -o ./ ${fastqc_report}
        """
}
