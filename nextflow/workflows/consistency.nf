#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


workflow CONSISTENCY {

    /* step 1. download data */
    GET_10X_DATA_FOR_CONSISTENCY_BENCHMARK(
        Channel.fromList(params.consistency.input_data.dataset),
        Channel.fromList(params.consistency.input_data.seed)
    )


}