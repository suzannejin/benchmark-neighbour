library(SingleCellExperiment)

pa <- argparser::arg_parser("Download 10X datasets from GEO dataset")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The GSE id of the dataset") 
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = 'the seed used for the simulation')
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "working directory")
pa <- argparser::add_argument(pa, "--download_helper", type = "character", help = "Download helper functions")
pa <- argparser::add_argument(pa, "--output", type = "character", help = "Path of the output file")

pa <- argparser::parse_args(pa)

set.seed(pa$seed)

source(pa$download_helper)

data_folder = pa$working_dir

# get data
message("Get ", pa$data_id)
sce <- data_loaders[[pa$data_id]]()
UMI <- as.matrix(assay(sce))

# downsampling
colsums <- colSums2(UMI)
downsample_proportion <- 5000 / median(colsums)
downsampled_UMI = as.matrix(scuttle::downsampleMatrix(UMI, prop = downsample_proportion, bycol = FALSE))

# keep only the expressed cells and genes
expressed_cells <- matrixStats::colSums2(downsampled_UMI) > 0
expressed_genes <- matrixStats::rowSums2(downsampled_UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells]
downsampled_UMI <- downsampled_UMI[expressed_genes, expressed_cells]

# save data
saveRDS(list(full = UMI, reduced = downsampled_UMI), pa$output)
