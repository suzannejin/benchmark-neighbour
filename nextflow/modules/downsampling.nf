
process GET_DATA_FOR_DOWNSAMPLING_BENCHMARK {

    input:
    each data_id
    each seed

    output:
    tuple val(data_id),
          val(seed),
          path("${data_id}_seed${seed}.rds"),
          emit: data

    script:
    def download_helper = params.download_helper
    """
    Rscript ${projectDir}/bin/downsampling_benchmark/download_deeply_sequenced_datasets.R \
        --data_id $data_id \
        --seed $seed \
        --working_dir . \
        --download_helper $download_helper
    """

    stub:
    def download_helper = params.download_helper
    """
    echo Rscript ${projectDir}/bin/downsampling_benchmark/download_deeply_sequenced_datasets.R \
        --data_id $data_id \
        --seed $seed \
        --working_dir . \
        --download_helper $download_helper
    touch ${data_id}_seed${seed}.rds
    """
}

process PLOT_SMARTSEQ3_DATA {

    input:
    tuple val(data_id),
          val(seed),
          path(data, stageAs: "data/*")

    output:
    tuple val(data_id),
          val(seed),
          path("${data_id}_seed${seed}.tsv")

    script:
    def transformation_helper = params.transformation_helper
    """
    mkdir cluster
    Rscript ${projectDir}/bin/downsampling_benchmark/plot_smartseq3_data.R \
        --data_id $data_id \
        --seed $seed \
        --working_dir . \
        --transformation_helper $transformation_helper
    """

    stub:
    def transformation_helper = params.transformation_helper
    """
    echo Rscript ${projectDir}/bin/downsampling_benchmark/plot_smartseq3_data.R \
        --data_id $data_id \
        --seed $seed \
        --working_dir . \
        --transformation_helper $transformation_helper
    mkdir cluster
    touch cluster/${data_id}_seed${seed}.tsv
    """
}

process TRANSFORM_DEEPLY_SEQUENCED_DATA {

    input:

    output:

    script:
}