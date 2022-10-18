library(rhdf5)
library(fastde)
library(tictoc)
library(future)

nworkers=64
future::plan("multicore", workers = nworkers, .cleanup=TRUE)

start_time <- Sys.time()
message("READ Sys.time start", start_time)
tic("READ brain 1.3M")
# f2 <- "~/data/SingleCell/1M_neurons_filtered_gene_bc_matrices_h5.h5"
f2 <- "~/scgc/data/1M_neurons_filtered_gene_bc_matrices_h5.h5"

sobject2 <- fastde::Read10X_h5_big(f2)
#sobject2
toc()
end_time <- Sys.time()
message("READ Sys.time end ", end_time, "duration", end_time - start_time)

nrows <- sobject2@Dim[1]
str(sobject2)

nclusters = 30

tic("GEN clusters")
# generate fake class labels.
clusters = 1:nclusters
labels2 <- sample(clusters, nrows, replace = TRUE)
toc()

start_time <- Sys.time()
message("FC SEXP Sys.time start", start_time)
tic("sparse fastde 64 FC SEXP")
# time and run BioQC
fastdefc3 <- fastde::ComputeFoldChangeSparse64SEXP(fastde::sp_transpose(sobject2), labels2, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = TRUE, threads = as.integer(nworkers))
toc()
end_time <- Sys.time()
message("FC SEXP Sys.time end ", end_time, "duration", end_time - start_time)


start_time <- Sys.time()
message("FC Sys.time start", start_time)
tic("sparse fastde 64 FC")
# time and run BioQC
fastdefc3 <- fastde::ComputeFoldChangeSparse64(fastde::sp_transpose(sobject2), labels2, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = TRUE, threads = as.integer(nworkers))
toc()
end_time <- Sys.time()
message("FC Sys.time end ", end_time, "duration", end_time - start_time)


colidx <- which(! colnames(fastdefc3) %in% c("cluster", "gene", "pct.1", "pct.2"))[1] 
start_time <- Sys.time()
message("FC Filter Sys.time start", start_time)
tic("sparse fastde 64 FC filter")
fc_mask <- fastde::FilterFoldChange(
    fastdefc3[[colidx]], fastdefc3$pct.1, fastdefc3$pct.2,
    init_mask = NULL,
    min_pct = 0.1, min_diff_pct = -Inf,
    logfc_threshold = 0.25, only_pos = FALSE,
    not_count = TRUE,
    threads = future::nbrOfWorkers())
toc()
end_time <- Sys.time()
message("FC filter Sys.time end ", end_time, "duration", end_time - start_time)




str(fastdefc3)

# run wilcox
tic("wilcox 1.3M")
wilcox_de2 <- fastde::FastFindAllMarkers64(sobject2, idents.clusters = labels2, test.use = 'fastwmw')
toc()
str(wilcox_de2)

# run ttest
#tic("ttest 1.3M")
#wilcox_de2 <- fastde::FastFindAllMarkers64(sobject2, idents.clusters = labels2, test.use = 'fast_t')
#toc()
#str(wilcox_de2)

message("FINISHED 4")
