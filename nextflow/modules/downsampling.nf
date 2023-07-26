
process GET_DATA_FOR_DOWNSAMPLING_BENCHMARK {

    input:
    each data_id
    each seed

    output:
    tuple val(data_id),
          val(seed),
          path("data/${data_id}.${seed}.rds"),
          emit: data

    script:
    def download_helper = params.download_helper
    """
    mkdir data/
    Rscript ${projectDir}/bin/downsampling_benchmark/download_deeply_sequenced_datasets.R \
        --data_id $data_id \
        --seed $seed \
        --working_dir . \
        --download_helper $download_helper \
        --output data/${data_id}.${seed}.rds 
    """

    stub:
    def download_helper = params.download_helper
    """
    echo Rscript ${projectDir}/bin/downsampling_benchmark/download_deeply_sequenced_datasets.R \
        --data_id $data_id \
        --seed $seed \
        --working_dir . \
        --download_helper $download_helper
    mkdir data/
    touch data/${data_id}.${seed}.rds
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
          path("cluster/${data_id}.${seed}.tsv")

    script:
    def transformation_helper = params.transformation_helper
    """
    mkdir cluster
    Rscript ${projectDir}/bin/downsampling_benchmark/plot_smartseq3_data.R \
        --data_id $data_id \
        --seed $seed \
        --working_dir . \
        --transformation_helper $transformation_helper \
        --input $data \
        --output cluster/${data_id}.${seed}.tsv
    """

    stub:
    def transformation_helper = params.transformation_helper
    """
    echo Rscript ${projectDir}/bin/downsampling_benchmark/plot_smartseq3_data.R \
        --data_id $data_id \
        --seed $seed \
        --working_dir . \
        --transformation_helper $transformation_helper \
        --input $data \
        --output cluster/${data_id}.${seed}.tsv
    mkdir cluster
    touch cluster/${data_id}.${seed}.tsv
    """
}

process TRANSFORM_DEEPLY_SEQUENCED_DATA {

    input:
    tuple val (data_id),
          val (seed),
          path (input_file, stageAs: 'data/*.rds'),
          val (transformation),
          val (alpha),
          val (knn),
          val (pca_dim),
          val (mode)
        
    output:
    tuple val (data_id),
          val (seed),
          val (transformation),
          val (knn),
          val (pca_dim),
          val (alpha),
          val (mode),
          path ("knn/${data_id}.${seed}.${transformation}.${knn}.${pca_dim}.${alpha}.${mode}.rds"),
          emit: knn
    tuple val (data_id),
          val (seed),
          val (transformation),
          val (knn),
          val (pca_dim),
          val (alpha),
          val (mode),
          path ("duration/${data_id}.${seed}.${transformation}.${knn}.${pca_dim}.${alpha}.${mode}.txt"),
          emit: duration

    script:
    def transformation_helper = params.transformation_helper

    """
    mkdir knn
    mkdir duration
    Rscript ${projectDir}/bin/downsampling_benchmark/transform_downsampling_data.R \
        --transformation $transformation \
        --knn $knn \
        --pca_dim $pca_dim \
        --alpha $alpha \
        --working_dir . \
        --transformation_helper $transformation_helper \
        --data_mode $mode \
        --input $input_file \
        --output_knn knn/${data_id}.${seed}.${transformation}.${knn}.${pca_dim}.${alpha}.${mode}.rds  \
        --output_duration duration/${data_id}.${seed}.${transformation}.${knn}.${pca_dim}.${alpha}.${mode}.txt
    """

    stub:
    def transformation_helper = params.transformation_helper

    """
    echo Rscript ${projectDir}/bin/downsampling_benchmark/transform_downsampling_data.R \
        --transformation $transformation \
        --knn $knn \
        --pca_dim $pca_dim \
        --alpha $alpha \
        --working_dir . \
        --transformation_helper $transformation_helper \
        --data_mode $mode \
        --input $input_file \
        --output_knn knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+${mode}.rds  \
        --output_duration duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+${mode}.txt
    mkdir knn
    mkdir duration
    touch knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+${mode}.rds
    touch duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+${mode}.txt
    """
}

process CALCULATE_DOWNSAMPLING_AGREEMENT {
    input:
    tuple val (data_id),
          val (seed),
          val (transformation),
          val (knn),
          val (pca_dim),
          val (alpha),
          val (mode)

    output:

    script:
    """
    mkdir results/downsampling/calculate
    Rscript ${projectDir}/bin/downsampling_benchmark/calculate_downsampling_agreement.R \
        --data_id $data_id \
        --seed $seed \
        --transformation $transformation \
        --knn $knn \
        --pca_dim $pca_dim \
        --alpha $alpha \
        --working_dir . \
        --full_knn_result_ids `ls results/downsampling/knn/ | grep full` \
        --reduced_knn_result_ids `ls results/downsampling/knn/ | grep reduced` \
        --result_path results/downsampling/calculate/${data_id}.${seed}.${transformation}.${knn}.${pca_dim}.${alpha}.${mode}.rds
    """
    stub:
    """
    mkdir results/downsampling/calculate
    echo Rscript ${projectDir}/bin/downsampling_benchmark/calculate_downsampling_agreement.R \
        --data_id $data_id \
        --seed $seed \
        --transformation $transformation \
        --knn $knn \
        --pca_dim $pca_dim \
        --alpha $alpha \
        --working_dir . \
        --full_knn_result_ids `ls results/downsampling/knn/ | grep full` \
        --reduced_knn_result_ids `ls results/downsampling/knn/ | grep reduced` \
        --result_path results/downsampling/calculate/${data_id}.${seed}.${transformation}.${knn}.${pca_dim}.${alpha}.${mode}.rds
        touch results/downsampling/calculate/stub.rds
    """

}