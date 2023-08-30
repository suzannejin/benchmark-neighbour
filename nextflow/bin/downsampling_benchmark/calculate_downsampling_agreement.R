library(tidyverse)

pa <- argparser::arg_parser("Calcualte the agreement between knns of the downsampled data and the standard of truth defined by the common knns identified by all transformations for the full data")
pa <- argparser::add_argument(pa, "--full_knns", type = "character", nargs = Inf, help = "Knn output file(s) for the full data") 
pa <- argparser::add_argument(pa, "--reduced_knns", type = "character", nargs = Inf, help = "Knn output file(s) for the reduced data")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "[Just for documentation purposes] A readable identifier of the data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "[Just for documentation purposes] The seed used to tame randomness")
pa <- argparser::add_argument(pa, "--transformations", type = "character", nargs = Inf, help = "[Just for documentation purposes] A readable identifier of the transformation(s)")
pa <- argparser::add_argument(pa, "--alphas", type = "character", nargs = Inf, help = "[Just for documentation purposes] Specification of the overdispersion(s)")
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "[Just for documentation purposes] The number of nearest neighbors considered")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "[Just for documentation purposes] The number of PCA dimensions used before KNN graph construction")
pa <- argparser::add_argument(pa, "--output", type = "character", help = "Output filename")
pa <- argparser::add_argument(pa, "--output_common_knns", type = "character", help = "Output filename to save the common knns")
pa <- argparser::parse_args(pa)


print(pa)

# check the same number of elements are given
stopifnot(length(pa$full_knns) == length(pa$reduced_knns))
stopifnot(length(pa$transformations) == length(pa$full_knns))
stopifnot(length(pa$alphas) == length(pa$full_knns))


# load knns
full_KNNs <- lapply(pa$full_knns, readRDS)
reduced_KNNs <- lapply(pa$reduced_knns, readRDS)
# Filter out negative controls
full_KNNs <- full_KNNs[! pa$transformations %in% c("raw_counts", "scaled_raw_counts")]

# Check that all knns have the same dimensions
stopifnot(nrow(full_KNNs[[1]]) == sapply(full_KNNs, nrow))
stopifnot(ncol(full_KNNs[[1]]) == sapply(full_KNNs, ncol))
stopifnot(nrow(full_KNNs[[1]]) == sapply(reduced_KNNs, nrow))
stopifnot(ncol(full_KNNs[[1]]) == sapply(reduced_KNNs, ncol))


# Define a set of reliable nearest neighbours
# as the set of knns of a cell that are common to all transformations
# on the deeply sequenced data (excluding negative controls)
n_cells <- nrow(full_KNNs[[1]])
common_knns <- lapply(seq_len(n_cells), function(idx){
  merged_nn <- lapply(full_KNNs, function(knn) knn[idx, ])
  purrr::reduce(merged_nn, intersect)
})
n_common_knns <- lengths(common_knns)
median_length_common_knns <- median(n_common_knns[n_common_knns > 1])


# Calculate the fraction of overlapping knns per cell
# between the ones found in the reduced set vs the standard of truth defined above
overlapping_knns_per_trans <- vapply(reduced_KNNs, function(knn){
  overlap <- sapply(seq_len(n_cells), function(idx){
    com_knn <- common_knns[[idx]]
    sum(knn[idx, ] %in% com_knn) 
  })
  mean(overlap[n_common_knns > 1])
}, numeric(1L))


# save results
res <- tibble(overlap = overlapping_knns_per_trans,
              median_length_common_knns = median_length_common_knns,
              dataset = pa$data_id, 
              seed = pa$seed, 
              transformation = pa$transformations, 
              alpha = pa$alphas,
              knn = pa$knn, 
              pca_dim = pa$pca_dim)
write_tsv(res, pa$output)