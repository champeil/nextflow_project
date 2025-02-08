# how to use?
    # nextflow run rnaseq.nf -c RNAseq.config -with-report report.html -with-trace -with-timeline timeline.html -with-dag dag.html
    # input: 
        # workflow .nf file
        # the corresponding config file
    # if use conda env, then need to add '-profile conda'
    # -c config is needed

# resume jobs
    # if failed while running, then use '-resume' parameters to resume the finished and paused jobs

# about the samplelist.csv
    # contain three column with header: sample_id, read1, read2
    # if single-end, then set read2 as "false"
    # the read1/read2 format: [sample_id]_[1/2]_fq.gz, in the absolute path

# standard of each process
    # input: input_dir, input_suffix, output_suffix, single_end, sample_id, path(input_file)
    # output: single_end, sample_id, path(output)

# about construct pipeline
    # construct pipeline is used to construct the pipeline automatically
        # according to the construct_info file
        # read up the construct_info file, and use the corresponding function, then echo to the pipeline file

# when place main.nf and nextflow.config together
    # then nextflow run main.nf will automatically load the nextflow.config in the same depth