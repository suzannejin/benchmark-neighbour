
process GET_DATA_FOR_DOWNSAMPLING_BENCHMARK {

    input:
    each data_id
    each seed

    output:
    tuple val(data_id),
          val(seed),
          path("data/${data_id}+${seed}.rds"),
          emit: data

    script:
    def download_helper = params.download_helper
    """
    mkdir data/
    Rscript ${projectDir}/bin/downsampling_benchmark/download_deeply_sequenced_datasets.R \
        --data_id $data_id \
        --seed $seed \
        --outdir data \
        --output data/${data_id}+${seed}.rds \
        --download_helper $download_helper 
    """

    stub:
    def download_helper = params.download_helper
    """
    echo Rscript ${projectDir}/bin/downsampling_benchmark/download_deeply_sequenced_datasets.R \
        --data_id $data_id \
        --seed $seed \
        --outdir data \
        --output data/${data_id}+${seed}.rds \
        --download_helper $download_helper 
    mkdir data/
    touch data/${data_id}+${seed}.rds
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
          path("cluster/${data_id}+${seed}.tsv")

    script:
    def transformation_helper = params.transformation_helper
    """
    mkdir cluster
    Rscript ${projectDir}/bin/downsampling_benchmark/plot_smartseq3_data.R \
        --data_id $data_id \
        --seed $seed \
        --input $data \
        --output cluster/${data_id}+${seed}.tsv \
        --transformation_helper $transformation_helper 
    """

    stub:
    def transformation_helper = params.transformation_helper
    """
    echo ${projectDir}/bin/downsampling_benchmark/plot_smartseq3_data.R \
        --data_id $data_id \
        --seed $seed \
        --input $data \
        --output cluster/${data_id}+${seed}.tsv \
        --transformation_helper $transformation_helper 
    mkdir cluster
    touch cluster/${data_id}+${seed}.tsv
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
          val (pca_dim)
        
    output:
    tuple val (data_id),
          val (seed),
          val (transformation),
          val (alpha),
          val (knn),
          val (pca_dim),
          path ("knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+full.rds"),
          path ("knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+reduced.rds"),
          emit: knn
    tuple val (data_id),
          val (seed),
          val (transformation),
          val (alpha),
          val (knn),
          val (pca_dim),
          path ("duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+full.txt"),
          path ("duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+reduced.txt"),
          emit: duration

    script:
    def transformation_helper = params.transformation_helper
    """
    mkdir knn
    mkdir duration
    Rscript ${projectDir}/bin/downsampling_benchmark/transform_downsampling_data.R \
        --input $input_file \
        --transformation $transformation \
        --alpha $alpha \
        --knn $knn \
        --pca_dim $pca_dim \
        --data_mode full \
        --output_knn knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+full.rds  \
        --output_duration duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+full.txt \
        --transformation_helper $transformation_helper 
    Rscript ${projectDir}/bin/downsampling_benchmark/transform_downsampling_data.R \
        --input $input_file \
        --transformation $transformation \
        --alpha $alpha \
        --knn $knn \
        --pca_dim $pca_dim \
        --data_mode reduced \
        --output_knn knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+reduced.rds  \
        --output_duration duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+reduced.txt \
        --transformation_helper $transformation_helper 
    """

    stub:
    def transformation_helper = params.transformation_helper
    """
    echo Rscript ${projectDir}/bin/downsampling_benchmark/transform_downsampling_data.R \
        --input $input_file \
        --transformation $transformation \
        --alpha $alpha \
        --knn $knn \
        --pca_dim $pca_dim \
        --data_mode full \
        --output_knn knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+full.rds  \
        --output_duration duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+full.txt \
        --transformation_helper $transformation_helper 
    echo Rscript ${projectDir}/bin/downsampling_benchmark/transform_downsampling_data.R \
        --input $input_file \
        --transformation $transformation \
        --alpha $alpha \
        --knn $knn \
        --pca_dim $pca_dim \
        --data_mode reduced \
        --output_knn knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+reduced.rds  \
        --output_duration duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+reduced.txt \
        --transformation_helper $transformation_helper 
    mkdir knn
    mkdir duration
    touch knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+full.rds 
    touch knn/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+reduced.rds
    touch duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+full.txt 
    touch duration/${data_id}+${seed}+${transformation}+${knn}+${pca_dim}+${alpha}+reduced.txt
    """
}

process CALCULATE_DOWNSAMPLING_AGREEMENT {

    input:
    tuple val (data_id),
          val (seed),
          val (transformations),
          val (alphas),
          val (knn),
          val (pca_dim),
          path(full_knns, stageAs: 'full_knns/*.rds'),
          path(reduced_knns, stageAs: 'reduced_knns/*.rds')

    output:
    tuple val (data_id),
          val (seed),
          val (knn),
          val (pca_dim),
          path("metric/${data_id}+${seed}+${knn}+${pca_dim}.tsv")

    script:
    def transformations2 = transformations.join(" ")
    def alphas2 = alphas.join(" ")
    """
    mkdir metric
    Rscript ${projectDir}/bin/downsampling_benchmark/calculate_downsampling_agreement.R \
        --full_knns $full_knns \
        --reduced_knns $reduced_knns \
        --data_id $data_id \
        --seed $seed \
        --transformations $transformations2 \
        --alphas $alphas2 \
        --knn $knn \
        --pca_dim $pca_dim \
        --output metric/${data_id}+${seed}+${knn}+${pca_dim}.tsv
    """

    stub:
    def transformations2 = transformations.join(" ")
    def alphas2 = alphas.join(" ")
    """
    mkdir metric
    echo Rscript ${projectDir}/bin/downsampling_benchmark/calculate_downsampling_agreement.R \
        --full_knns $full_knns \
        --reduced_knns $reduced_knns \
        --data_id $data_id \
        --seed $seed \
        --transformations $transformations2 \
        --alphas $alphas2 \
        --knn $knn \
        --pca_dim $pca_dim \
        --output metric/${data_id}+${seed}+${knn}+${pca_dim}.tsv
    touch metric/${data_id}+${seed}+${knn}+${pca_dim}.tsv
    """

}