# FastDE

This is an R package for performing differential gene expression analysis.   It is designed as a drop-in replacement for the Seurat FindAllMarkers function.

The goal of the library is to improve computational performance and large scale data handling.

Currently, Seurat supports 9 algorithms.   FastDE supports student's t-test and Wilcox Rank Sum test.


## Installation

### From source.
TO BE COMPLETED.

	BiocManager::install(c("argparser", "SCnorm", "sva", "limma", "SingleCellExperiment", "edgeR"), lib=.libPaths()[1])
	
	BiocManager::install(c("future", "future.apply", "stringi", "reshape2", "magrittr", "tidyr", "fastcluster", "diceR"), lib=.libPaths()[1])
	
	BiocManager::install(c("proftools", "profvis", "tictoc"), lib=.libPaths()[1])
