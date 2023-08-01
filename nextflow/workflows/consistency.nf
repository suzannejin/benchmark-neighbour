#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GET_10X_DATA_FOR_CONSISTENCY_BENCHMARK; PLOT_10X_DATA; TRANSFORM_CONSISTENCY_DATA; CALCULATE_10X_CONSISTENCY } from "../modules/consistency.nf"


workflow CONSISTENCY {


    /* step 1. download data */
    GET_10X_DATA_FOR_CONSISTENCY_BENCHMARK( Channel.fromList(params.consistency.input_data.dataset) )

    // plot t-SNE and UMAP
    PLOT_10X_DATA( GET_10X_DATA_FOR_CONSISTENCY_BENCHMARK.out.data )


    /* step 2. transform consistency data and generate a KNN graph */
    GET_10X_DATA_FOR_CONSISTENCY_BENCHMARK.out.data
        .combine( Channel.fromList(params.consistency.input_data.seed) )
        .combine( Channel.fromList(params.consistency.knn_construction.transformations) )
        .map{ it -> [
            it[0],
            it[1],
            it[2],
            it[3].name, 
            it[3].alpha   // TODO check alpha parameters
        ]}
        .transpose(by:4)
        .combine( Channel.fromList(params.consistency.knn_construction.knn) )
        .combine( Channel.fromList(params.consistency.knn_construction.pca) )
        .set{ data2transform2knn }  
    TRANSFORM_CONSISTENCY_DATA(data2transform2knn)


    /* step 3. calculate consistency metric */
    CALCULATE_10X_CONSISTENCY( TRANSFORM_CONSISTENCY_DATA.out.knn )
    
}