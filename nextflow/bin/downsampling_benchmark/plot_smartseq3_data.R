
pa <- argparser::arg_parser("Take simulated data and generate the tSNE and UMAP for plotting")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The input data id")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "The seed that was used for the downsampling")
pa <- argparser::add_argument(pa, "--input", type = "character", help = "The input data path")
pa <- argparser::add_argument(pa, "--output", type = "character", help = "The output data path")
pa <- argparser::add_argument(pa, "--transformation_helper", type = "character", help = "Transformation helper functions")
pa <- argparser::parse_args(pa)

print(pa)

source(pa$transformation_helper)

# Get ground truth and simulated data
counts_full <- readRDS(pa$input)$full
counts_reduced <- readRDS(pa$input)$reduced

# get average cell totals
sf_full <- MatrixGenerics::colSums2(counts_full)
sf_full <- sf_full / mean(sf_full)
sf_reduced <- MatrixGenerics::colSums2(counts_reduced)
sf_reduced <- sf_reduced / mean(sf_reduced)

# log expression
log_counts_full <- logp1_fnc(counts_full, sf_full, alpha = FALSE)
log_counts_reduced <- logp1_fnc(counts_reduced, sf_reduced, alpha = FALSE)

# run pca
pca_log_counts_full <-  BiocSingular::runPCA(t(log_counts_full), rank = 2, get.rotation = FALSE, BSPARAM = BiocSingular::FastAutoParam())$x
pca_log_counts_reduced <-  BiocSingular::runPCA(t(log_counts_reduced), rank = 2, get.rotation = FALSE, BSPARAM = BiocSingular::FastAutoParam())$x

# run tsne
tsne_log_counts_full <- scater::calculateTSNE(log_counts_full)
tsne_log_counts_reduced <- scater::calculateTSNE(log_counts_reduced)

# do clustering
clustering <- bluster::clusterRows(t(log_counts_full), bluster::KNNGraphParam(), full = TRUE)
if(length(unique(clustering$clusters)) > 15){
  clusters <- bluster::mergeCommunities(clustering$objects$graph, clustering$clusters, number = 15)
}else{
  clusters <- clustering$clusters
}

# save data
res <- data.frame(name = pa$data_id, seed=pa$seed, cluster = clusters, 
                  col_sums_full = MatrixGenerics::colSums2(counts_full),
                  col_sums_reduced = MatrixGenerics::colSums2(counts_reduced),
                  tsne_log_counts_full_axis1 = tsne_log_counts_full[,1],       tsne_log_counts_full_axis2 = tsne_log_counts_full[,2],
                  tsne_log_counts_reduced_axis1 = tsne_log_counts_reduced[,1], tsne_log_counts_reduced_axis2 = tsne_log_counts_reduced[,2],
                  pca_log_counts_full_axis1 = pca_log_counts_full[,1],         pca_log_counts_full_axis2 = pca_log_counts_full[,2],
                  pca_log_counts_reduced_axis1 = pca_log_counts_reduced[,1],   pca_log_counts_reduced_axis2 = pca_log_counts_reduced[,2])
write.table(res, pa$output, sep = "\t", row.names = FALSE, quote = FALSE)