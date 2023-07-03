#!/usr/bin/env nextflow

GET_10X_DATA_FOR_CONSISTENCY_BENCHMARK {

    input:
    each val(dataset)
    each val(seed)

    output:
    

    script:
    def download_helper = params.download_helper
    """
    Rscript ${projectDir}/bin/consistency_benchmark/download_10X_dataset.R \
        --dataset $dataset \
        --data_folder . \
        --download_helper ${projectDir}/bin/consistency_benchmark/download_helper.R
    """

    stub:
    def download_helper = params.download_helper
    """
    echo Rscript ${projectDir}/bin/consistency_benchmark/download_10X_dataset.R \
        --dataset $dataset \
        --data_folder . \
        --download_helper ${projectDir}/bin/consistency_benchmark/download_helper.R
    """
}