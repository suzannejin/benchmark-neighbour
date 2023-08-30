#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GET_DATA_FOR_DOWNSAMPLING_BENCHMARK } from "../modules/downsampling.nf"
include { PLOT_SMARTSEQ3_DATA                 } from "../modules/downsampling.nf"
include { TRANSFORM_DEEPLY_SEQUENCED_DATA     } from "../modules/downsampling.nf"
include { CALCULATE_DOWNSAMPLING_AGREEMENT    } from "../modules/downsampling.nf"

workflow DOWNSAMPLING {

    /* step 1. download data */
    GET_DATA_FOR_DOWNSAMPLING_BENCHMARK(
        Channel.fromList(params.downsampling.input_data.dataset),
        Channel.fromList(params.downsampling.input_data.seed),
    )

    // plot pca and tsne
    PLOT_SMARTSEQ3_DATA( GET_DATA_FOR_DOWNSAMPLING_BENCHMARK.out.data )

    /* step 2. transform data and generate a KNN graph */
    GET_DATA_FOR_DOWNSAMPLING_BENCHMARK.out.data
        .combine( Channel.fromList(params.downsampling.knn_construction.transformations) )
        .map{ it -> [
            it[0],
            it[1],
            it[2],
            it[3].name, 
            it[3].alpha 
        ]}
        .transpose(by:4)
        .combine( Channel.fromList(params.downsampling.knn_construction.knn) )
        .combine( Channel.fromList(params.downsampling.knn_construction.pca) )
        .set{ data2transform2knn }
    TRANSFORM_DEEPLY_SEQUENCED_DATA(data2transform2knn)

    /* step 3. calculate downsampling agreement */
    TRANSFORM_DEEPLY_SEQUENCED_DATA.out.knn
        .groupTuple(by:[0,1,4,5])
        .set{knn2overlap}
    CALCULATE_DOWNSAMPLING_AGREEMENT(knn2overlap)
}