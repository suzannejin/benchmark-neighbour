#We want a unique id, with nf we don't need different names for data and result
pa <- argparser::arg_parser("Take a matrix and generate a KNN graph")
pa <- argparser::add_argument(pa, "--transformation", type = "character", help = "The name of the transformation")
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "The number of k nearest neighbors that are considered")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "The dimensions for the pca transformation")
pa <- argparser::add_argument(pa, "--alpha", type = "character", default = "FALSE", help = "The alpha parameter. Ignored by some transformations.")
pa <- argparser::add_argument(pa, "--transformation_helper", type = "character", help = "Transformation helper functions")
pa <- argparser::add_argument(pa, "--data_mode", type = "character", help = "Either 'full' or 'reduced'")
pa <- argparser::add_argument(pa, "--input", type = "character", help = "Path of the input file.")
pa <- argparser::add_argument(pa, "--output_duration", type = "character", help = "Path of the output file (.txt)")
pa <- argparser::add_argument(pa, "--output_knn", type = "character", help = "Path of the knn output file (.rds)")


pa <- argparser::parse_args(pa)

print(pa)
stopifnot(pa$data_mode %in% c("full", "reduced"))


# get a list called `all_transformations` that contains 
# all transformations as ready to call functions
# Furthermore, it initializes the `make_knn_graph` function
source(pa$transformation_helper)


######### Start Transformation #########
#Read UMI matrix

if(pa$data_mode == "full"){
  UMI <- readRDS(pa$input)$full
}else{
  UMI <- readRDS(pa$input)$reduced
}
expressed_cells <- matrixStats::colSums2(UMI) > 0
expressed_genes <- matrixStats::rowSums2(UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells]

#Parsing alpha
alpha <- pa$alpha
if(pa$alpha == "global"){
  alpha <- "global"
}else if(! is.na(suppressWarnings(readr::parse_double(pa$alpha, na = character(0L))))){
  alpha <- readr::parse_double(pa$alpha)
}else if(! is.na(suppressWarnings(readr::parse_logical(pa$alpha, na = character(0L))))){
  alpha <- readr::parse_logical(pa$alpha)
}else{
  stop("Cannot parse alpha=", alpha)
}

sf <- MatrixGenerics::colSums2(UMI)
sf <- sf / mean(sf)

#Transform data and make KNN graph
duration <- system.time({
  trans_dat <- all_transformations[[pa$transformation]](UMI, sf, alpha)
  KNN <- make_knn_graph(pa$transformation, trans_dat, pa$pca_dim, pa$knn)
})

write.table(data.frame(name = names(duration), seconds = as.vector(duration)),
            pa$output_duration,
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(KNN, pa$output_knn)
