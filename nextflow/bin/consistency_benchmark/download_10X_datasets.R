library(SingleCellExperiment)

pa <- argparser::arg_parser("Download 10X datasets from GEO dataset")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The GSE id of the dataset") 
pa <- argparser::add_argument(pa, "--outdir", type = "character", help = "Folder where downloaded data should be stored")
pa <- argparser::add_argument(pa, "--output", type = "character", help = "Output filename") 
pa <- argparser::add_argument(pa, "--download_helper", type = "character", help = "Download helper functions")
pa <- argparser::parse_args(pa)

source(pa$download_helper)
data_folder <- pa$outdir

message("Get ", pa$dataset)
sce <- data_loaders[[pa$data_id]]()
UMI <- as.matrix(assay(sce))

message("Saving ", pa$data_id)
saveRDS(UMI, pa$output)
