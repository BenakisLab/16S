#!/usr/bin/env Rscript
# Author: Adam Sorbie
# Date: 18/04/24
# Version 1.1.0

if(!require("pacman")){
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
pacman::p_load(dada2, optparse, parallel, ggplot2, tictoc, plotly,htmlwidgets, R.utils)


option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path of read files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="dada2_qc", 
              help="output folder", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Functions

RIGHT <- function(x,n){
  substring(x,nchar(x)-n+1)
}

# Start of script

# Convert paths to absolute
opt$path <- getAbsolutePath(opt$path)
opt$out <- getAbsolutePath(opt$out)

# check all paths have trailing forward slash 
if (RIGHT(opt$path, 1) != "/") {
  opt$path <- paste0(opt$path, "/")
}
if (RIGHT(opt$out, 1) != "/") {
  opt$out <- paste0(opt$out, "/")
}


if (is.null(opt$path)){
  print_help(opt_parser)
  stop("Please provide path", 
       call.=FALSE)
}

#check path exists
if (!dir.exists(opt$path)){
  stop("Supplied directory does not exist", call. = FALSE)
}

# time analysis
tic()
print(paste("QC STARTING", Sys.time(), sep=" "))

path <- opt$path 
dir.create(opt$out, showWarnings = F)
setwd(opt$path)


fnFs <- sort(list.files(pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pattern="_R2_001.fastq.gz", full.names = TRUE))

# get sample names from forward reads 
sample.names <- sapply(strsplit(basename(fnFs), "_S1_L001_R1"), `[`, 1)



print("Plotting quality profiles")

qc_F <- plotQualityProfile(fnFs, aggregate=TRUE) + 
  ggtitle("Forward")
int_F <- ggplotly(qc_F)

if (length(fnRs)!=0){
  qc_R <- plotQualityProfile(fnRs, aggregate=TRUE) + 
    ggtitle("Reverse")
  int_R <- ggplotly(qc_R)
  
  plot_out <- subplot(int_F, int_R) %>% 
    layout(title = "Forward                Reverse")
  htmlwidgets::saveWidget(widget = plot_out, file=paste(opt$out,"quality_profiles.html", sep="/"))
} 




print(paste("QC COMPLETED", Sys.time(), sep=" "))

sink(paste(opt$out, "parameter_log.txt", sep="/"))
print(paste0("Filepath: ", opt$path))
print(paste0("Output: ", opt$out))
sessionInfo()
toc()
closeAllConnections() 
