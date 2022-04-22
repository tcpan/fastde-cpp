unlink(".RData")

library("devtools")
#library("proftools")

options(width=120)

# https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

print(thisFile())
devtools::load_all(file.path(dirname(thisFile()), "..", "..", "R"))


tic("TOTAL RUN TIME")
p <- build_workflow_argparser()
# Parse the command line arguments

in_args <- commandArgs(trailingOnly = TRUE)
print(in_args)
argv <- argparser::parse_args(p, in_args)

validate_workflow_args(p, argv)
print_workflow_args(argv)

workflow_function <- if (argv$integ == "beer") {
    beer_cluster
} else if (argv$integ == "harmony") {
    harmony_cluster
} else if (argv$integ == "seurat") {
    seurat_cluster
} else if (argv$integ == "combat") {
    combat_cluster
} else {
    undefined_cluster
}

#prof_name=paste(argv$out_dir, paste("rprof", argv$integ, argv$qc, argv$vis, argv$clust_algorithm, argv$gen_markers, "out", sep="."), sep="/")
#Rprof(filename=prof_name, line.profiling=TRUE, gc.profiling = TRUE)

if (argv$workers > 1) {
    options(future.globals.maxSize = g_maxsize_limit)
    print(paste("available cores orig", future::availableCores(), sep=' '))
    options(mc.cores = argv$workers)
    print(paste("set cores to", argv$workers, sep=' '))
    print(paste("available cores new", future::availableCores(), sep=' '))
    future::plan("multicore", workers = argv$workers, .cleanup=TRUE)
} else {
    future::plan("sequential")
}

workflow_function(argv$root_dir, argv$data_file_csv,
            argv$out_dir, argv$qc, argv$vis, argv$img, argv$test,
            argv$save, argv$inc, argv$exc, argv$clust_algorithm,
            argv$load, argv$gen_markers, argv$dot_markers,
            argv$subset_list_file)

#Rprof()
toc()

## Explicitly close multisession workers by switching plan
if (argv$workers > 1) {
    future::plan("sequential")  
}
#pd <- readProfileData(prof_name)
#funSum <- funSummary(pd)
#funSum <- funSum[order(funSum$self.pct, decreasing=TRUE), ]
#print(funSum[1:20, ])
#callSum <- callSummary(pd)
#callSum <- callSum[order(callSum$self.pct, decreasing=TRUE), ]
#print(callSum[1:30, ])
#srcSum <- srcSummary(pd)
#srcSum <- srcSum[order(srcSum$total.pct, decreasing=TRUE), ]
#print(srcSummary(pd)[1:30, ])
#hotPaths(pd, total.pct=1.0)
# summaryRprof(filename=prof_name, lines="both", diff=TRUE)

