#!/usr/bin/env nextflow

process GET_10X_DATA_FOR_CONSISTENCY_BENCHMARK {

    input:
    each dataset

    output:
    tuple val(dataset), 
          path("${dataset}.rds"),
          emit: data
    tuple val(dataset), 
          path("GSM*/*"),
          emit: gsm

    script:
    def download_helper = params.download_helper
    """
    Rscript ${projectDir}/bin/consistency_benchmark/download_10X_datasets.R \
        --dataset $dataset \
        --data_folder . \
        --download_helper $download_helper
    """

    stub:
    def download_helper = params.download_helper
    """
    echo Rscript ${projectDir}/bin/consistency_benchmark/download_10X_datasets.R \
        --dataset $dataset \
        --data_folder . \
        --download_helper $download_helper
    touch ${dataset}.rds
    mkdir GSM
    touch GSM/test
    """
}

process PLOT_10X_DATA {

    input:
    tuple val(data_id),
          path(dataset, stageAs: "data/*")

    output:
    tuple val(data_id), 
          path("cluster/${data_id}.tsv")

    script:
    def transformation_helper = params.transformation_helper
    """
    mkdir cluster
    Rscript ${projectDir}/bin/consistency_benchmark/plot_10x_data.R \
        --data_id $data_id \
        --working_dir . \
        --transformation_helper $transformation_helper
    """

    stub:
    def transformation_helper = params.transformation_helper
    """
    echo Rscript ${projectDir}/bin/consistency_benchmark/plot_10x_data.R \
        --data_id $data_id \
        --working_dir . \
        --transformation_helper $transformation_helper
    mkdir cluster
    touch cluster/${data_id}.tsv
    """
}

process TRANSFORM_CONSISTENCY_DATA {

    input: 
    tuple val(data_id),
          path(dataset, stageAs: 'data/*'),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim),
          val(seed)
    
    output:
    tuple val(data_id),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim),
          val(seed),
          path("knn/${data_id}.rds"),
          emit: knn
    tuple val(data_id),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim),
          val(seed),
          path("duration/${data_id}.txt"),
          emit: duration

    script:
    def transformation_helper = params.transformation_helper
    """
    mkdir duration
    mkdir knn
    Rscript ${projectDir}/bin/consistency_benchmark/transform_consistency_data.R \
        --transformation $transformation \
        --data_id $data_id \
        --knn $knn \
        --pca_dim $pca_dim \
        --alpha $alpha \
        --seed $seed \
        --transformation_helper $transformation_helper \
        --working_dir .
    """

    stub:
    def transformation_helper = params.transformation_helper
    """
    echo Rscript ${projectDir}/bin/consistency_benchmark/transform_consistency_data.R \
        --transformation $transformation \
        --data_id $data_id \
        --knn $knn \
        --pca_dim $pca_dim \
        --alpha $alpha \
        --seed $seed \
        --transformation_helper $transformation_helper \
        --working_dir .
    mkdir knn
    mkdir duration
    touch knn/${data_id}.rds
    touch duration/${data_id}.txt
    """
}

process CALCULATE_10X_CONSISTENCY {

    input:
    tuple val(data_id),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim),
          val(seed),
          path(data, stageAs: "knn/*")

    output:
    tuple val(data_id),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim),
          val(seed),
          path("metric/${data_id}.tsv")

    script:
    """
    mkdir metric
    Rscript ${projectDir}/bin/consistency_benchmark/calculate_10X_consistency.R \
        --transformation $transformation \
        --data_id $data_id \
        --knn $knn \
        --pca_dim $pca_dim \
        --alpha $alpha \
        --seed $seed \
        --working_dir .
    """

    stub:
    """
    echo Rscript ${projectDir}/bin/consistency_benchmark/calculate_10X_consistency.R \
        --transformation $transformation \
        --data_id $data_id \
        --knn $knn \
        --pca_dim $pca_dim \
        --alpha $alpha \
        --seed $seed \
        --working_dir .
    mkdir metric
    touch metric/${data_id}.tsv
    """
}