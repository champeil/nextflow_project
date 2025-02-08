#!/usr/bin/env nextflow
// This process is for featurecount quantification
// author: laojp
// time: 2025.02.06
// position: SYSUCC bioinformatic platform
// usage
//      input: input_dir, input_suffix, output_suffix, sample_id, map_bam

nextflow.enable.dsl=2

process process_featurecount {
    label "featurecount"
    tag "${sample_id}"
    publishDir "${params.result_dir}/featurecount", mode: 'copy'

    input:
        tuple val(input_dir), val(input_suffix), val(output_suffix), val(single_end), val(sample_id), path(map_bam)
    output:
        val("${sample_id}")
        val("${single_end}")
        path("${sample_id}${output_suffix}.txt"), emit: featurecount_output
        path("${sample_id}${output_suffix}.txt.summary")

    script:
        def featurecount_cmd = ""
        if (single_end) {
            // single  end
            featurecount_cmd = """
                ${params.software_featurecount}/featureCounts \
                    -T $task.cpus \
                    -t exon \
                    -g gene_id \
                    -a ${params.reference_gtf} \
                    -o ${sample_id}${output_suffix}.txt \
                    ${map_bam}
            """
        } else {
            // paired end
            featurecount_cmd = """
                ${params.software_featurecount}/featureCounts \
                    -p --countReadPairs \
                    -T $task.cpus \
                    -t exon \
                    -g gene_id \
                    -a ${params.reference_gtf} \
                    -o ${sample_id}${output_suffix}.txt \
                    ${map_bam}
            """
        }

        // run featurecount
        """
            #!/usr/bin/env bash
            ${featurecount_cmd}
        """
}

process process_featurecount_merge {
    label "featurecount_merge"
    tag "featurecount_merge"
    publishDir "${params.result_dir}/featurecount", mode: 'copy'

    input:
        tuple val(input_dir), val(input_suffix), val(output_suffix), path(featurecount_res)
    output:
        path("featurecount_merge.txt")

    script:
        """
            #!/usr/bin/env bash
            cat ${featurecount_res[0]} | grep -v "#" | awk -v FS="\t" '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6}' > featurecount_merge_temp.txt
            for feature_count in ${featurecount_res}; do
                name=\$(basename \${feature_count} ${input_suffix}.txt)
                paste featurecount_merge_temp.txt <(grep -v "#" \${feature_count} | awk -v name=\${name} -v FS="\t" '{if(NR==1){\$7=name};print \$7}') > temp ;
                mv temp featurecount_merge_temp.txt
            done
            mv featurecount_merge_temp.txt featurecount_merge.txt
        """
}