#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

if (params.run_consistency_benchmark){
    include { CONSISTENCY } from "./workflows/consistency"
}
if (params.run_downsampling_benchmark){
    include { DOWNSAMPLING } from "./workflows/downsampling.nf"
}

workflow {
    if (params.run_consistency_benchmark){
        CONSISTENCY()
    }
    if (params.run_downsampling_benchmark){
        DOWNSAMPLING()
    }
}