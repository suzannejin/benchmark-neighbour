#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GET_DATA_FOR_DOWNSAMPLING_BENCHMARK } from "../modules/downsampling.nf"


workflow DOWNSAMPLING {

    /* step 1. download data */
    GET_DATA_FOR_DOWNSAMPLING_BENCHMARK(
        Channel.fromList(params.downsampling.input_data.dataset),
        Channel.fromList(params.downsampling.input_data.seed),
    )

    // plot pca and tsne
    PLOT_SMARTSEQ3_DATA( GET_DATA_FOR_DOWNSAMPLING_BENCHMARK.out.data )

    

}