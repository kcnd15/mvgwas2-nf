#!/usr/bin/env Rscript

library(MTAR)

library(optparse)

# command line arguments
option_list <- list(
  make_option(c("-d", "--datadir"), type = "character",
              help = "data directory [default '%default']", 
              default = "../../data/mtar", metavar = "FILE"),
  make_option(c("-c", "--covariates"), type = "character",
              help = "Covariate file (indId, [covariates]), [default '%default']", 
              default = "cov.csv", metavar = "FILE"),
  make_option(c("-t", "--traits"), type = "character",
              help = "Traits file (indId, [traits]), [default '%default']", 
              default = "traits.csv", metavar = "FILE"),
  make_option(c("-g", "--genotypes"), type = "character",
              help = "Genotypes file (indId, [genotype]), [default '%default']", 
              default = "geno.csv", metavar = "FILE"),
  make_option(c("-v", "--verbose"), action = "store_true", 
              help = "[default %default]", 
              default = TRUE),
  make_option(c("-s", "--sample"), action = "store_true", 
              help = "[default %default]", 
              default = FALSE)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Step 1: Prepare Summary Statistics for Multiple Traits

cat("MTAR - prepare summary statistics\n\n")

# 1. Given the Individual-Level Data

if (opt$sample) {
  
  data("rawdata")
  attach(rawdata)
  
  obs.stat <- Get_UV_from_data(traits = traits.dat,
                               covariates = cov.dat,
                               genotype = geno.dat,
                               covariance = TRUE)
  
  print(obs.stat$U)
  
  detach(rawdata)
  stop()
}

## 2. Input files

traits.f <- paste0(opt$datadir, "/", opt$traits)
cov.f <- paste0(opt$datadir, "/", opt$covariates) 
geno.f <- paste0(opt$datadir, "/", opt$genotypes) 

if ( is.null(traits.f) || is.null(geno.f) || is.null(cov.f)) {
  print_help(opt_parser)
  stop("Missing/not found I/O files")
}

cat("traits     :", traits.f, "\n")
cat("covariates :", cov.f, "\n")
cat("genotypes  :", geno.f, "\n\n")

# read input files and create the input structure for MTAR

# read csv, convert dataframe to list of matrix
get_list_matrix <- function(csv_path) {
  
  # read csv-file with header into dataframe
  study.csv = read.csv2(file = csv_path , header = TRUE, sep=",", dec=".")
  
  # convert to matrix
  study.matrix = as.matrix(study.csv)
  
  # set rownames
  study.matrix.rows = nrow(study.matrix)
  row_names = paste0("SUBJ", 1:study.matrix.rows)
  dimnames(study.matrix)[[1]] = row_names
  study.dat = list(Study1=study.matrix)
  
  # return list of matrix
  return (study.dat)
}

# input file
# cov.f = "/home/kcan/UOC/tfm/mvgwas2-nf/data/mtar/cov.csv"
# traits.f = "/home/kcan/UOC/tfm/mvgwas2-nf/data/mtar/traits.csv"
# geno.f = "/home/kcan/UOC/tfm/mvgwas2-nf/data/mtar/geno.csv"

# convert input data to list of matrix
cov.dat = get_list_matrix(cov.f)
traits.dat = get_list_matrix(traits.f)
geno.dat = get_list_matrix(geno.f)

# call MTAR function
obs.stat <- MTAR::Get_UV_from_data(traits = traits.dat,
                             covariates = cov.dat,
                             genotype = geno.dat,
                             covariance = TRUE)

# show results
obs.stat$U
