#!/usr/bin/env nextflow
// This process is for star mapping
// author: laojp
// time: 2025.02.06
// position: SYSUCC bioinformatic platform
// usage
//      input: input_dir, input_suffix, output_suffix, sample_id, reads

nextflow.enable.dsl=2

process process_star {
    label "STAR"
    tag "$sample_id"
    publishDir "${params.result_dir}/star", mode: 'copy'

    input:
        tuple val(input_dir), val(input_suffix), val(output_suffix), val(single_end), val(sample_id), path(reads)
    output:
        tuple val("${sample_id}"),
            val("${single_end}"),
            path("${sample_id}${output_suffix}*Aligned.sortedByCoord.out.bam"),
            path("${sample_id}${output_suffix}*Aligned.toTranscriptome.out.bam"),
            path("${sample_id}${output_suffix}*out.tab"),
            path("${sample_id}${output_suffix}*out"), emit: star_output
    
    script:
        def star_cmd = ""
        if (single_end) {
            // single end
            star_cmd = """
                ${params.software_star}/STAR \
                    --runThreadN ${task.cpus} \
                    --genomeDir ${params.reference_star} \
                    --readFilesIn ${reads[0]} \
                    --readFilesCommand zcat \
                    --outFileNamePrefix ${sample_id}${output_suffix} \
                    --outSAMtype BAM SortedByCoordinate \
                    --quantMode TranscriptomeSAM GeneCounts \
                    --outReadsUnmapped Fastx/
            """
        } else {
            // paired end
            star_cmd = """
                ${params.software_star}/STAR \
                    --runThreadN ${task.cpus} \
                    --genomeDir ${params.reference_star} \
                    --readFilesIn ${reads[0]} ${reads[1]} \
                    --readFilesCommand zcat \
                    --outFileNamePrefix ${sample_id}${output_suffix} \
                    --outSAMtype BAM SortedByCoordinate \
                    --quantMode TranscriptomeSAM GeneCounts \
                    --outReadsUnmapped Fastx/
            """
        }

        // run star
        """
            #!/usr/bin/env bash
            ${star_cmd}
        """


}
