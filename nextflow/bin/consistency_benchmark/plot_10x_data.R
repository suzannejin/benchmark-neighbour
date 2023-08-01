
pa <- argparser::arg_parser("Take simulated data and generate the tSNE and UMAP for plotting")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The name of the dataset")
pa <- argparser::add_argument(pa, "--input", type = "character", help = "Path to input data (format .rds)")
pa <- argparser::add_argument(pa, "--output", type = "character", help = "Path to output file")
pa <- argparser::add_argument(pa, "--transformation_helper", type = "character", help = "Transformation helper functions")
pa <- argparser::parse_args(pa)
print(pa)

source(pa$transformation_helper)

# read count data
counts <- readRDS(pa$input)

# get average cell totals
sf <- MatrixGenerics::colSums2(counts)
sf <- sf / mean(sf)

# get log counts
log_counts <- logp1_fnc(counts, sf, alpha = FALSE)

# run pca
pca_log_counts <- BiocSingular::runPCA(t(log_counts), rank = 20, get.rotation = FALSE, BSPARAM = BiocSingular::FastAutoParam())$x

# run tsne
tsne_log_counts <- scater::calculateTSNE(log_counts)

# run clustering
clustering <- bluster::clusterRows(pca_log_counts, bluster::KNNGraphParam(), full = TRUE)
if(length(unique(clustering$clusters)) > 15){
  clusters <- bluster::mergeCommunities(clustering$objects$graph, clustering$clusters, number = 15)
}else{
  clusters <- clustering$clusters
}

# save data
res <- data.frame(name = pa$data_id, cluster = clusters, col_sums = MatrixGenerics::colSums2(counts),
                  tsne_log_counts_axis1 = tsne_log_counts[,1], tsne_log_counts_axis2 = tsne_log_counts[,2],
                  pca_log_counts_axis1 = pca_log_counts[,1],   pca_log_counts_axis2 = pca_log_counts[,2])
write.table(res, pa$output, sep = "\t", row.names = FALSE, quote = FALSE)