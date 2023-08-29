library(tidyverse)

pa <- argparser::arg_parser("Take a matrix and generate a KNN graph")
pa <- argparser::add_argument(pa, "--input", type = "character", help = "Path to input knn .rds file") 
pa <- argparser::add_argument(pa, "--output", type = "character", help = "Path to output file") 
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "[Just for documentation purposes] The seed used to tame randomness")
pa <- argparser::add_argument(pa, "--transformation", type = "character", help = "[Just for documentation purposes] A readable identifier of the transformation")
pa <- argparser::add_argument(pa, "--alpha", type = "character", help = "[Just for documentation purposes] Specification of the overdispersion.")
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "[Just for documentation purposes] The number of nearest neighbors considered")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "[Just for documentation purposes] The number of PCA dimensions used before KNN graph construction")
pa <- argparser::parse_args(pa)

print(pa)


KNNs <- readRDS(pa$input)
stopifnot(all(dim(KNNs[[1]]) == dim(KNNs[[2]])))
n_genes <- nrow(KNNs[[1]])


cons <- mean(sapply(seq_len(n_genes), function(gene_idx){
  length(intersect(KNNs[[1]][gene_idx,], KNNs[[2]][gene_idx,]))
}))

res <- tibble(mean_overlap = cons, dataset = pa$data_id, seed = pa$seed,
              pca_dim = pa$pca_dim, knn = pa$knn, transformation = pa$transformation, alpha = pa$alpha)

write_tsv(res, pa$output)