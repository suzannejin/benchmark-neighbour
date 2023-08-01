#!/usr/bin/env nextflow

process GET_10X_DATA_FOR_CONSISTENCY_BENCHMARK {

    input:
    each data_id

    output:
    tuple val(data_id), 
          path("data/${data_id}.rds"),
          emit: data
    tuple val(data_id), 
          path("data/GSM*/*"),
          emit: gsm

    script:
    def download_helper = params.consistency_download_helper
    """
    mkdir data
    Rscript ${projectDir}/bin/consistency_benchmark/download_10X_datasets.R \
        --data_id $data_id \
        --outdir data \
        --output data/${data_id}.rds \
        --download_helper $download_helper
    """

    stub:
    def download_helper = params.consistency_download_helper
    """
    mkdir data
    echo Rscript ${projectDir}/bin/consistency_benchmark/download_10X_datasets.R \
        --data_id $data_id \
        --outdir data \
        --output data/${data_id}.rds \
        --download_helper $download_helper
    touch data/${data_id}.rds
    mkdir data/GSM
    touch data/GSM/test
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
        --input $dataset \
        --output cluster/${data_id}.tsv \
        --transformation_helper $transformation_helper
    """

    stub:
    def transformation_helper = params.transformation_helper
    """
    mkdir cluster
    echo Rscript ${projectDir}/bin/consistency_benchmark/plot_10x_data.R \
        --data_id $data_id \
        --input $dataset \
        --output cluster/${data_id}.tsv \
        --transformation_helper $transformation_helper
    touch cluster/${data_id}.tsv
    """
}

process TRANSFORM_CONSISTENCY_DATA {

    input: 
    tuple val(data_id),
          path(dataset, stageAs: 'data/*'),
          val(seed),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim)
    
    output:
    tuple val(data_id),
          val(seed),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim),
          path("knn/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.rds"),
          emit: knn
    tuple val(data_id),
          val(seed),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim),
          path("duration/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.txt"),
          emit: duration

    script:
    def transformation_helper = params.transformation_helper
    """
    mkdir duration
    mkdir knn
    Rscript ${projectDir}/bin/consistency_benchmark/transform_consistency_data.R \
        --input $dataset \
        --output_duration duration/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.txt \
        --output_knn knn/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.rds \
        --seed $seed \
        --transformation $transformation \
        --alpha $alpha \
        --knn $knn \
        --pca_dim $pca_dim \
        --transformation_helper $transformation_helper 
    """

    stub:
    def transformation_helper = params.transformation_helper
    """
    mkdir duration
    mkdir knn
    echo Rscript ${projectDir}/bin/consistency_benchmark/transform_consistency_data.R \
        --input $dataset \
        --output_duration duration/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.txt \
        --output_knn knn/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.rds \
        --seed $seed \
        --transformation $transformation \
        --alpha $alpha \
        --knn $knn \
        --pca_dim $pca_dim \
        --transformation_helper $transformation_helper 
    touch duration/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.txt
    touch knn/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.rds
    """
}

process CALCULATE_10X_CONSISTENCY {

    input:
    tuple val(data_id),
          val(seed),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim),
          path(data, stageAs: "knn/*")

    output:
    tuple val(data_id),
          val(seed),
          val(transformation),
          val(alpha),
          val(knn),
          val(pca_dim),
          path("metric/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.tsv")

    script:
    """
    mkdir metric
    Rscript ${projectDir}/bin/consistency_benchmark/calculate_10X_consistency.R \
        --input $data \
        --output metric/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.tsv \
        --seed $seed \
        --transformation $transformation \
        --alpha $alpha \
        --knn $knn \
        --pca_dim $pca_dim 
    """

    stub:
    """
    mkdir metric
    echo Rscript ${projectDir}/bin/consistency_benchmark/calculate_10X_consistency.R \
        --input $data \
        --output metric/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.tsv \
        --seed $seed \
        --transformation $transformation \
        --alpha $alpha \
        --knn $knn \
        --pca_dim $pca_dim 
    touch metric/${data_id}+${seed}+${transformation}+${alpha}+${knn}+${pca_dim}.tsv
    """
}