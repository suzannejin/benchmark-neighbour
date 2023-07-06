#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

if (params.run_consistency_benchmark){
    include { CONSISTENCY } from "./workflows/consistency"
}

workflow {
    if (params.run_consistency_benchmark){
        CONSISTENCY()
    }
}