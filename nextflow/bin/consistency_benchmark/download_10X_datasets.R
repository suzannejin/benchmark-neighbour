library(SingleCellExperiment)

pa <- argparser::arg_parser("Download 10X datasets from GEO dataset")
pa <- argparser::add_argument(pa, "--dataset", type = "character", help = "The GSE id of the dataset") 
pa <- argparser::add_argument(pa, "--data_folder", type = "character", help = "Folder where data should be stored")
pa <- argparser::add_argument(pa, "--download_helper", type = "character", help = "Download helper functions")
pa <- argparser::parse_args(pa)

source(pa$download_helper)
data_folder <- pa$data_folder

message("Get ", pa$dataset)
sce <- data_loaders[[pa$dataset]]()
UMI <- as.matrix(assay(sce))

message("Saving ", pa$dataset)
saveRDS(UMI, file.path(data_folder, paste0(pa$dataset, '.rds')))
