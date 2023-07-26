#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GET_DATA_FOR_DOWNSAMPLING_BENCHMARK } from "../modules/downsampling.nf"
include { PLOT_SMARTSEQ3_DATA } from "../modules/downsampling.nf"
include { TRANSFORM_DEEPLY_SEQUENCED_DATA } from "../modules/downsampling.nf"

workflow DOWNSAMPLING {

    /* step 1. download data */
    GET_DATA_FOR_DOWNSAMPLING_BENCHMARK(
        Channel.fromList(params.downsampling.input_data.dataset),
        Channel.fromList(params.downsampling.input_data.seed),
    )

    // plot pca and tsne
    PLOT_SMARTSEQ3_DATA( GET_DATA_FOR_DOWNSAMPLING_BENCHMARK.out.data )

    /* step 2. transform consistency data and generate a KNN graph */
    GET_DATA_FOR_DOWNSAMPLING_BENCHMARK.out.data
        .combine( Channel.fromList(params.consistency.knn_construction.transformations) )
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
        .combine( Channel.fromList(["full", "reduced"]) )
        .set{ data2transform2knn }
    TRANSFORM_DEEPLY_SEQUENCED_DATA(data2transform2knn)

    // calculate
    CALCULATE_DOWNSAMPLING_AGREEMENT(TRANSFORM_DEEPLY_SEQUENCED_DATA.out.knn)
}