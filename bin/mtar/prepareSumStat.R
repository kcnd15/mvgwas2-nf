#!/usr/bin/env Rscript

library(MTAR)

library(optparse)

# command line arguments
option_list <- list(
  make_option(c("-c", "--covariates"), type = "character",
              help = "Covariate file (indId, [covariates])", 
              metavar = "FILE"),
  make_option(c("-t", "--traits"), type = "character",
              help = "Covariate file (indId, [covariates])", 
              metavar = "FILE"),
  make_option(c("-g", "--genotypes"), type = "character",
              help = "Covariate file (indId, [covariates])", 
              metavar = "FILE"),
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

cat("MTAR - prepare summary statistics\n")

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

traits.f <- opt$traits 
cov.f <- opt$covariates
geno.f <- opt$genotypes

if ( is.null(traits.f) || is.null(geno.f) || is.null(cov.f)) {
  print_help(opt_parser)
  stop("Missing/not found I/O files")
}

cat("covariates:", opt$covariates, "\n")
