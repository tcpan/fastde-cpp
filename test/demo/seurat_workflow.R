#' @importFrom magrittr "%>%"
library(tictoc)

create_seurat_object <- function(data_combmat, data_batch_names,
                     normalize = TRUE, project = "INTGSC", dim_name = "dataset",
                     mincells = 5, nfeat = 2000, npcs = 20) {
    data_sobj <- Seurat::CreateSeuratObject(counts = data_combmat,
                               project = project, min.cells = mincells)
    data_sobj@meta.data[dim_name] <- data_batch_names
    if (normalize == TRUE) {
        data_sobj <- data_sobj %>% Seurat::NormalizeData(verbose = FALSE)
    }
    data_sobj <- data_sobj %>% Seurat::ScaleData(verbose = FALSE)
    data_sobj <- data_sobj %>% Seurat::FindVariableFeatures(
        selection.method = "vst", nfeatures = nfeat, verbose = FALSE)
    #
    if (npcs > 0) {
        data_sobj <- data_sobj %>% Seurat::RunPCA(
            pc.genes = (data_sobj)@var.genes, npcs = npcs, verbose = FALSE)
    }
    data_sobj
}


create_seurat_objects_list <- function(mat10x_list, data_names,
                           short_names, dim_name = "dataset", normalize = TRUE,
                           min_cells = 3,  min_feats = 200, nfeats = 2000) {

    du_lapply(seq_len(length(mat10x_list)), function(i) {
        cmatx <- mat10x_list[[i]]
        if (is.na(colnames(cmatx)) ||
            (length(colnames(cmatx)) != dim(cmatx)[2])) {
            stop(paste(short_names[i], "matrix has bad  colnames",
                 length(colnames(cmatx)), " ",
            dim(cmatx)[2], "\n"))
        }
        # colnames(cmatx) <- stringr::str_c(short_names[i], colnames(cmatx),
        #                                  sep = "-")
        sobj <- Seurat::CreateSeuratObject(counts = cmatx,
                    project = data_names[i], min.cells = min_cells,
                    min.features = min_feats)
        if (normalize) {
            sobj <- Seurat::NormalizeData(sobj)
        }
        sobj <- Seurat::FindVariableFeatures(sobj, selection.method = "vst",
                    nfeatures = nfeats)
        sobj@meta.data[dim_name] <- c(rep(short_names[i], ncol(cmatx)))
        sobj
    })
}

load_10x_seurat_objs <- function(base_dir, dir_paths, data_names,
                        short_names, dim_name = "dataset", min_cells = 3,
                        min_feats = 200, nfeats = 2000) {

    mat10x_list <- load_10x_matrices(base_dir, dir_paths, data_names)
    mat10x_sobjects <- create_seurat_objects_list(
                        mat10x_list, data_names, short_names, dim_name, TRUE,
                        min_cells, min_feats, nfeats)
    cat("Seurat Objects loaded for all datasets\n")
    mat10x_sobjects
}

load_10x_seurat_objs_union <- function(base_dir, dir_paths, data_names,
                        short_names, dim_name = "dataset", min_cells = 3,
                        min_feats = 200, nfeats = 2000) {

    mat10x_list <- load_10x_matrices_union(base_dir, dir_paths, data_names)
    mat10x_sobjects <- create_seurat_objects_list(
                        mat10x_list, data_names, short_names, dim_name, TRUE,
                        min_cells, min_feats, nfeats)
    cat("Seurat Objects loaded for all datasets\n")
    mat10x_sobjects
}

qcload_10x_seurat_objs <- function(base_dir, dir_paths, data_names,
          short_names, qc_function, inc_list_file, exc_list_file, qc_option,
          subset_list = NULL,
          dim_name = "dataset", nfeats = 2000, min_cells = 2, min_feats = 200) {

    mat10x_list <- qcload_10x_matrices(base_dir, dir_paths, data_names,
                                qc_function, inc_list_file, exc_list_file,
                                qc_option)
    if (!is.null(subset_list)) {
        nmat_list <- length(mat10x_list)
        tmp_mat_list <- du_lapply(seq_len(nmat_list),
            function(i) {
                cnames <- intersect(colnames(mat10x_list[[i]]),
                                             subset_list)
                if (length(cnames) > 0) {
                    cat("Subset for ", short_names[i], ": ", 
                        dim(mat10x_list[[i]]), length(subset_list),
                        length(cnames),
                        dim(mat10x_list[[i]][, cnames]), "\n")
                    mat10x_list[[i]][, cnames]
                } else {
                    cat("Subset for ", short_names[i], ": ",
                        length(subset_list), length(cnames),
                        dim(mat10x_list[[i]]), "\n")
                    mat10x_list[[i]]
                }
            })
        mat10x_list <- tmp_mat_list
    }

    # Should Normalize ? sobj = NormalizeData(sobj)
    mat10x_sobjects <- create_seurat_objects_list(
                        mat10x_list, data_names, short_names, dim_name, FALSE,
                        min_cells, min_feats, nfeats)
    cat("Seurat Objects loaded for all datasets\n")
    mat10x_sobjects
}

qcload_10x_seurat_objs_union <- function(base_dir, dir_paths, data_names,
     short_names, qc_function, inc_list_file, exc_list_file, qc_option,
     subset_list=NULL,
     dim_name = "dataset", nfeats = 2000, min_cells = 2, min_feats = 200) {

    mat10x_list <- qcload_10x_matrices_union(base_dir, dir_paths, data_names,
                                qc_function, inc_list_file, exc_list_file,
                                qc_option)
    if (!is.null(subset_list)) {
        # TODO: subset list
        nmat_list <- length(mat10x_list)
        mat10x_list <- du_lapply(seq_len(nmat_list),
            function(i) {
                cnames <- intersect(colnames(mat10x_list[i]),
                                             subset_list)
                if (length(cnames) > 0)
                    mat10x_list[[i]][, cnames]
                else
                    mat10x_list[[i]]
            })
    }

    # Should Normalize ? sobj = NormalizeData(sobj)
    mat10x_sobjects <- create_seurat_objects_list(
                        mat10x_list, data_names, short_names, dim_name, FALSE,
                        min_cells, min_feats, nfeats)
    cat("Seurat Objects loaded for all datasets\n")
    mat10x_sobjects
}

combined_seurat_object <- function(data_mlist, short_names, normalize = TRUE,
                                project = "ATHSC", dim_name = "dataset",
                                mincells = 5, nfeat = 2000, npcs = 20) {
    clst <- combined_expr_matrix(data_mlist, short_names)
    data_combmat <- clst[[1]]
    data_batch_names <- clst[[2]]
    create_seurat_object(data_combmat, data_batch_names, normalize, project,
        dim_name, mincells, nfeat, npcs)
}


integrate_seurat_objects <- function(seurat_obj_lst, anc_feat, ndims = 1:30) {
    integ_anchors <- Seurat::FindIntegrationAnchors(
        object.list = seurat_obj_lst, anchor.features = anc_feat, dims = ndims)
    Seurat::IntegrateData(anchorset = integ_anchors, dims = ndims)
}

#" @export
seurat_allqc_plot <- function(data_dir, out_dir, out_prefix, img_option) {
    scr_data <- Seurat::Read10X(data_dir = data_dir)
    scrj <- Seurat::CreateSeuratObject(counts = scr_data, project = "SCRNA",
                                       min.cells = 3, min.features = 200)
    scrj[["percent.mt"]] <- Seurat::PercentageFeatureSet(scrj,
                                                         pattern = "^ATMG")
    p1 <- Seurat::VlnPlot(scrj, features = c("nFeature_RNA", "nCount_RNA",
                                             "percent.mt"), ncol = 3)
    p2 <- Seurat::FeatureScatter(scrj, feature1 = "nCount_RNA",
                                 feature2 = "percent.mt")
    p3 <- Seurat::FeatureScatter(scrj, feature1 = "nCount_RNA",
                                 feature2 = "nFeature_RNA")
    p4 <- Seurat::CombinePlots(plots = list(p2, p3))
    qfname <- paste(out_dir, "/", out_prefix, "qcvln.", img_option, sep = "")
    ggplot2::ggsave(qfname, p1)
    ggplot2::ggsave(paste(out_dir, "/", out_prefix, " fscatter.", img_option,
                          sep = ""), p4, width = 14, height = 7)

    scrj <- Seurat::NormalizeData(scrj, normalization.method = "LogNormalize",
                                  scale.factor = 10000)
    scrj <- Seurat::FindVariableFeatures(scrj, selection.method = "vst",
                                         nfeatures = 2000)
    top10 <- head(Seurat::VariableFeatures(scrj), 10)
    p5 <- Seurat::VariableFeaturePlot(scrj)
    p6 <- Seurat::LabelPoints(plot = p5, points = top10, repel = TRUE)
    p7 <- Seurat::CombinePlots(plots = list(p5, p6))
    ggplot2::ggsave(paste(out_dir, "/", out_prefix, "vfeat.", img_option,
                          sep = ""), p7, width = 14, height = 7)
    scrj <- Seurat::ScaleData(scrj)
    paste("Data Scaled \n")
    scrj <- Seurat::RunPCA(scrj,
                           features = Seruat::VariableFeatures(object = scrj))
    p8 <- Seurat::DimPlot(scrj, reduction = "pca")
    ggplot2::ggsave(paste(out_dir, "/", out_prefix, "pca.", img_option,
                          sep = ""),
                    p8)
}

seurat_fscatter <- function(scrj, fname) {
    scrj[["percent.mcg"]] <- Seurat::PercentageFeatureSet(scrj,
                                                          pattern = "^AT[MC]G")
    p2 <- Seurat::FeatureScatter(scrj, feature1 = "nCount_RNA",
                                 feature2 = "percent.mcg")
    p3 <- Seurat::FeatureScatter(scrj, feature1 = "nCount_RNA",
                                 feature2 = "nFeature_RNA")
    p4 <- Seurat::CombinePlots(plots = list(p2, p3))
    ggplot2::ggsave(fname, p4, width = 14, height = 7)
}

seurat_ba_qcplot <- function(data_dir, all_thrs, dirx, out_dir, out_prefix,
                             image_option = "png") {
    scr_data <- Seurat::Read10X(data_dir = data_dir)
    scrj <- Seurat::CreateSeuratObject(counts = scr_data, project = data_dir)
    scrj[["percent.mcg"]] <- Seurat::PercentageFeatureSet(scrj,
                                                          pattern = "^AT[MC]G")
    fname <- paste(out_dir, dirx, paste(out_prefix, "-seurat-before-qc.",
                                        image_option, sep = ""), sep = "/")
    seurat_fscatter(scrj, fname)


    mcg_threshold <- all_thrs[1]
    min_cells <- all_thrs[2]
    min_features <- all_thrs[3]
    cat("MCG minc", data_dir, mcg_threshold, min_cells, "\n")
    scrj <- Seurat::CreateSeuratObject(counts = scr_data, project = data_dir,
                         min.cells = all_thrs[2], min.features = min_features)
    scrj[["percent.mcg"]] <- Seruat::PercentageFeatureSet(scrj,
                                                          pattern = "^AT[MC]G")
    # srch = subset(pbmc, subset = nCount_RNA > 200 & nFeature_RNA < 2500 &
    # percent.mcg < mcg_threshold)
    # scrj = subset(scrj, subset = nCount_RNA > 200 &
    #               percent.mcg < mcg_threshold)
    expr1 <- Seurat::FetchData(object = scrj, vars = "nFeature_RNA")
    scrj1 <- scrj[, expr1 > 200]
    expr2 <- Seurat::FetchData(object = scrj1, vars = "percent.mcg")
    scrj2 <- scrj1[, expr2 < mcg_threshold]

    fname <- paste(out_dir, dirx, paste(out_prefix, "-seuerat-after-qc.",
                                        image_option, sep = ""), sep = "/")
    seurat_fscatter(scrj2, fname)
}


seurat_dim_plot <- function(sobj, reduce_by, group = "dataset", split = NULL,
                     label = FALSE, width = 6, height = 4, out_file = NULL) {
    if (!is.null(out_file)) {
        options(repr.plot.height = height, repr.plot.width = width)
    }
    px <- Seurat::DimPlot(object = sobj, reduction = reduce_by, pt.size = 0.1,
                          group.by = group, split.by = split, label = label)
    if (!is.null(out_file)) {
        ggplot2::ggsave(out_file, px, width = width, height = height)
    }
    px
}

seurat_violin_plot <- function(sobj, feats, group = "dataset",
                           width = 6, height = 4,  out_file = NULL) {
    if (!is.null(out_file)) {
        options(repr.plot.height = height, repr.plot.width = width)
    }
    px <- Seurat::VlnPlot(object = sobj, features = feats, group.by = group,
                          pt.size = 0.1)
    if (!is.null(out_file)) {
        ggplot2::ggsave(out_file, px, width = width, height = height)
    }
    px
}

seurat_dim_violin_plot <- function(sobj, reduce_by, feats, group,
                                   out_file = NULL) {
    options(repr.plot.height = 5, repr.plot.width = 12)
    p1 <- seurat_dim_plot(sobj, reduce_by, group)
    p2 <- seurat_violin_plot(sobj, feats, group)
    p3 <- cowplot::plot_grid(p1, p2)
    if (!is.null(out_file)) {
        ggplot2::ggsave(out_file, p3, width = 12, height = 5)
    }
    p3
}

cluster_umap_seurat <- function(data_sobj, reduce_by, resolution = 0.5,
                                dims = 1:20, clust_algorithm = 1) {
    tic("[TIME] Seurat Run UMAP")
    data_sobj <- data_sobj %>% Seurat::RunUMAP(reduction = reduce_by,
                                               dims = dims)
    toc()
    tic("[TIME] Seurat Find Neighbors")
    data_sobj <- data_sobj %>% Seurat::FindNeighbors(reduction = reduce_by,
                                                    dims = dims)
    toc()
    tic("[TIME] Seurat Find Clusters")
    data_sobj <- data_sobj %>% Seurat::FindClusters(resolution = resolution,
                                                    algorithm = clust_algorithm)
    toc()
    data_sobj <- data_sobj %>% identity()
    data_sobj
}


cluster_tsne_seurat <- function(data_sobj, reduce_by, resolution = 0.5,
                                dims = 1:20, clust_algorithm = 1) {
    tic("[TIME] Seurat Run TSNE")
    data_sobj <- data_sobj %>% Seurat::RunTSNE(reduction = reduce_by,
                                               dims = dims)
    toc()
    tic("[TIME] Seurat Find Neighbors")
    data_sobj <- data_sobj %>% Seurat::FindNeighbors(reduction = reduce_by,
                                                    dims = dims)
    toc()
    tic("[TIME] Seurat Find Clusters")
    data_sobj <- data_sobj %>% Seurat::FindClusters(resolution = resolution,
                                                    algorithm = clust_algorithm)
    toc()
    data_sobj <- data_sobj %>% identity()
    data_sobj
}



seurat_combined_umap <- function(seurat_object, out_dir, reduce_by = "pca",
                        dims = 1:30, clust_algorithm = 1,
                        img_option = "png", out_suffix = "") {
    #
    seurat_object <- cluster_umap_seurat(seurat_object, reduce_by, 0.5, dims,
                                     clust_algorithm)
    #
    print(seurat_object@reductions)
    # seurat_dim_plot(seurat_object, reduce_by = "umap", group = "dataset",
    #          split = "dataset", width = 10, height = 4,
    #          out_file = paste(out_dir, paste(reduce_by,
    #             "-seurat-umap-grouped", out_suffix, ".",
    #                  img_option,  sep = ""), sep = "/"))
    #
    seurat_dim_plot(
        seurat_object, reduce_by = "umap", group = NULL, label = TRUE,
        width = 6, height = 4, out_file = paste(out_dir, paste(reduce_by,
            "-seurat-umap-integrated", out_suffix, ".",
                     img_option,  sep = ""), sep = "/"))
    seurat_object
}


seurat_combined_tsne <- function(seurat_object, out_dir, reduce_by = "pca",
                                 dims = 1:30, clust_algorithm = 1,
                                 img_option = "png", out_suffix = "") {
    #
    seurat_object <- cluster_tsne_seurat(seurat_object, reduce_by, 0.5, dims,
                                     clust_algorithm)
    #
    print(seurat_object@reductions)
    # seurat_dim_plot(seurat_object, reduce_by = "tsne", group = "dataset",
    #          split = "dataset", width = 10, height = 4,
    #          out_file = paste(out_dir, paste(reduce_by,
    #             "-seurat-tsne-grouped", out_suffix, ".",
    #                  img_option,  sep = ""), sep = "/"))
    #
    seurat_dim_plot(
        seurat_object, reduce_by = "tsne", group = NULL, label = TRUE,
        width = 6, height = 4, out_file = paste(out_dir, paste(reduce_by,
            "-seurat-tsne-integrated.", out_suffix, ".",
                     img_option,  sep = ""), sep = "/"))
    seurat_object
}

seurat_integrate <- function(root_dir, data_file, out_dir, qc_option,
                           vis_option, img_option, save_option, save_suffix,
                           inc_list_file, exec_list_file, clust_algorithm,
                           subset_list=NULL) {
    tic("[TIME] Seurat Load")

    data_df <- read.csv(data_file, header = TRUE, stringsAsFactors = FALSE)
    print(head(data_df))
    expt_dir_paths <- data_df[, "dir.paths"]
    short_names <- data_df[, "short.names"]
    project_names <- data_df[, "project.names"]

    seurat_obj_lst <- if (!is.na(qc_option)) {
        qcload_10x_seurat_objs_union(root_dir, expt_dir_paths,
                            project_names, short_names, qc_normalize_matrix,
                            inc_list_file, exec_list_file, qc_option,
                            subset_list)
    } else {
        load_10x_seurat_objs(root_dir, expt_dir_paths, project_names,
                                short_names, subset_list)
    }
    # for (s in seurat_obj_lst) {
    #     print(s@meta.data)
    # }
    toc()

    tic("[TIME] Seurat Integrate")
    if (length(seurat_obj_lst) > 1) {
        integrated_sobject <- integrate_seurat_objects(seurat_obj_lst, 16000)
        Seurat::DefaultAssay(integrated_sobject) <- "integrated"
    } else {
        integrated_sobject <- seurat_obj_lst[[1]]
    }
    print(integrated_sobject)
    toc()

    tic("[TIME] Seurat Scale")
    # print(integrated_sobject@meta.data)
    # Run the standard workflow for visualization and clustering
    integrated_sobject <- Seurat::ScaleData(integrated_sobject,
                                              verbose = FALSE)
    toc()

    tic("[TIME] Seurat PCA")
    integrated_sobject <- Seurat::RunPCA(integrated_sobject, npcs = 30,
                                           verbose = FALSE)
    toc()

    tic("[TIME] Seurat Visualize")
    if (!is.na(vis_option)) {
        if (vis_option == "umap") {
            integrated_sobject <- seurat_combined_umap(integrated_sobject,
                out_dir, "pca", 1:30, clust_algorithm, img_option, save_suffix)
        } else if (vis_option == "tsne") {
            integrated_sobject <- seurat_combined_tsne(integrated_sobject,
                out_dir, "pca", 1:30, clust_algorithm, img_option, save_suffix)
        }
    }
    toc()

    tic("[TIME] Seurat Save")
    if (save_option) {
        save(integrated_sobject, file = paste(out_dir, "/",
            "Seurat", save_suffix, ".Rdata", sep = ""))
    }
    toc()
    integrated_sobject
}

seurat_find_markers <- function(integrated_sobject, out_dir, vis_option,
                             img_option, test_option, out_prefix,
                             all_markers_file, dot_markers_file) {
    tic("[TIME] Seurat Find Markers")
    if (!(is.na(all_markers_file))) {
        cat("Generating All Markers ... \n")
        if ((test_option == "fastwmw") || (test_option == "bioqc") || (test_option == "fast_t")) {
            mkdf <- fastde::FastFindAllMarkers(integrated_sobject,
                                                test.use = test_option)
        } else {
            mkdf <- Seurat::FindAllMarkers(integrated_sobject,
                                        test.use = test_option)
        }
        write.table(mkdf, all_markers_file, row.names = FALSE, sep = "\t")
    }
    toc()

    tic("[TIME] Seurat Plot Markers")
    if (!is.na(dot_markers_file)) {
        mdf <- read.table(dot_markers_file, stringsAsFactors = F, sep = "\t")
        for (dxf in mdf$V1) {
            gdf <- read.table(dxf, header = T, sep = "\t", stringsAsFactors = F)
            genes_plot <- gdf$ID
            pfx <- gsub("\\.tsv", "", basename(dxf))
            npresent <- genes_plot %in% rownames(integrated_sobject)
            cat(dxf, pfx, sum(npresent), sum(!npresent), "\n")
            if (sum(npresent) > 0) {
                dot_fname <- paste(out_dir, paste(pfx, "-dot-markers-seurat-",
                            out_prefix, ".", img_option, sep = ""), sep = "/")
                px <- Seurat::DotPlot(integrated_sobject,
                                      features = genes_plot) +
                            ggplot2::ggtitle(pfx)
                ggplot2::ggsave(dot_fname, px, width = 16, height = 8)
            }
        }
    }
    toc()

}

#' @export
seurat_cluster <- function(root_dir, data_file, out_dir, qc_option, vis_option,
                           img_option, test_option, save_option,
                           inc_list_file, exec_list_file,
                           clust_algorithm, robject_file,
                           gen_markers, dot_markers_file,
                           subset_list_file=NA) {

    subset_list <- if (is.na(subset_list_file)) {
        NULL
    } else {
        read.table(subset_list_file, header = FALSE, 
                   stringsAsFactors = FALSE)$V1
    }

    save_suffix <- paste("-", qc_option, "-c", clust_algorithm, sep = "")
    integrated_sobject <- if (!is.na(robject_file)) {
        cat("Loading object file ", robject_file, "...\n")
        load(robject_file)
        integrated_sobject
    } else {
        seurat_integrate(root_dir, data_file, out_dir, qc_option, vis_option,
                         img_option, save_option, save_suffix,
                         inc_list_file, exec_list_file, clust_algorithm,
                         subset_list)
    }
    print(integrated_sobject@reductions)

    all_markers_file <- if (gen_markers) {
         paste(out_dir, paste(vis_option, "-seurat-markers-", save_suffix,
                              ".tsv", sep = ""),
               sep = "/")
    } else {
     NA
    }

    seurat_find_markers(integrated_sobject, out_dir, vis_option, img_option,
        test_option, save_suffix, all_markers_file, dot_markers_file)

    integrated_sobject
}


#' @export
seurat_custom_cluster <- function(root_dir, data_file, out_dir, cluster_file,
                  qc_option, vis_option, img_option, test_option, save_option,
                  inc_list_file, exec_list_file, clust_algorithm,
                  robject_file, gen_markers, dot_markers_file,
                  subset_list_file=NA) {

    subset_list <- if (is.na(subset_list_file)) {
        NULL
    } else {
        read.table(subset_list_file, header = FALSE, 
                   stringsAsFactors = FALSE)$V1
    }


    save_suffix <- paste("-custom-", qc_option, "c",
                          clust_algorithm, sep = "")
    integrated_sobject <- if (!is.na(robject_file)) {
        cat("Loading object file ", robject_file, "...\n")
        load(robject_file)
        integrated_sobject
    } else {
        seurat_integrate(root_dir, data_file, out_dir, qc_option, "none",
                       img_option, save_option, save_suffix,
                       inc_list_file, exec_list_file, clust_algorithm,
                       subset_list)
    }
    out_suffix <- ""
    clust_df <- read.table(cluster_file, header = FALSE,
                           stringsAsFactors = FALSE,
                           sep = "\t", col.names = c("Cell", "Cluster"))
    integrated_sobject$cust_grp <- clust_df$Cluster
    seurat_dim_plot(integrated_sobject, reduce_by = "umap", group = "cust_grp",
                    split = "dataset", width = 10, height = 4,
                    out_file = paste(out_dir,
                         paste("Custom-beer-custom-umap-grouped", out_suffix,
                               ".", img_option, sep = ""), sep = "/"))
    Seurat::Idents(integrated_sobject) <- "cust_grp"
    all_markers_file <- if (gen_markers) {
        paste(out_dir, paste(vis_option, "-beer-custom-markers.tsv", sep = ""),
              sep = "/")
    } else {
        NA
    }
    seurat_find_markers(integrated_sobject, out_dir, vis_option, img_option,
                      test_option, save_suffix,
                      all_markers_file, dot_markers_file)

    integrated_sobject

}
