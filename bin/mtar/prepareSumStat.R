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
              help = "Genotypes matrix file (indId, [genotype]), [default '%default']", 
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
get_list_matrix <- function(csv_path, sep=",") {
  
  # read csv-file with header into dataframe
  study.csv = read.csv2(file = csv_path , header = TRUE, sep=sep, dec=".")
  
  # convert to matrix
  nr = nrow(study.csv)
  nc = ncol(study.csv)
  
  # get row names = first column of csv
  study.csv.rownames.list = as.list(study.csv[1])
  study.csv.rownames = study.csv.rownames.list[[1]]
  study.csv.data = study.csv[-1] # remove first column which is the row name
  
  # convert strings to factors
  
  study.matrix = as.matrix(study.csv, nrow=nr, ncol=nc)
  
  rownames(study.matrix) = study.csv.rownames
  dn = dimnames(study.matrix)
  # dimnames(study.matrix)[[1]] = study.csv.rownames
  
  # set rownames
  study.matrix.rows = nrow(study.matrix)
  row_names = paste0("SUBJ", 1:study.matrix.rows)
  dimnames(study.matrix)[[1]] = row_names
  study2.dat = list(Study1=study.matrix)
  study.dat = list(Study1=study.csv)
  
  # return list of matrix
  return (study2.dat)
}

# read csv, convert dataframe to list of matrix
get_covariates <- function(csv_path, sep=",") {
  
  # read csv-file with header into dataframe
  study.csv = read.csv2(file = csv_path , header = TRUE, sep=sep, dec=".")
  
  # get row names = first column of csv
  study.csv.rownames.list = as.list(study.csv[1])
  study.csv.rownames = study.csv.rownames.list[[1]]
  study.csv.data = study.csv[-1] # remove first column which is the row name
  
  # convert strings to factors
  # study.sex_factor = as.factor(study.csv.data[2])
  # study.sex = as.integer(study.sex_factor)
  sv = as.vector(study.csv.data[,2])
  svi = as.integer(as.factor(sv))
  study.csv.data[,2] = svi
  
  # convert to matrix
  nr = nrow(study.csv.data)
  nc = ncol(study.csv.data)
  study.matrix = as.matrix(study.csv.data, nrow=nr, ncol=nc)
  
  # set rownames
  rownames(study.matrix) = study.csv.rownames
  
  # return list of matrix
  study.dat = list(Study1=study.matrix)
  return (study.dat)
}

# read csv, convert dataframe to list of matrix
get_traits <- function(csv_path, sep=",") {
  
  # read csv-file with header into dataframe
  study.csv = read.csv2(file = csv_path , header = TRUE, sep=sep, dec=".")
  
  # get row names = first column of csv
  study.csv.rownames.list = as.list(study.csv[1])
  study.csv.rownames = study.csv.rownames.list[[1]]
  study.csv.data = study.csv[-1] # remove first column which is the row name
  
  # convert strings to factors
  # study.sex_factor = as.factor(study.csv.data[2])
  # study.sex = as.integer(study.sex_factor)
  
  # not necessary for traits, no factors
  # sv = as.vector(study.csv.data[,2])
  # svi = as.integer(as.factor(sv))
  # study.csv.data[,2] = svi
  
  # convert to matrix
  nr = nrow(study.csv.data)
  nc = ncol(study.csv.data)
  study.matrix = as.matrix(study.csv.data, nrow=nr, ncol=nc)
  
  # set rownames
  rownames(study.matrix) = study.csv.rownames
  
  # return list of matrix
  study.dat = list(Study1=study.matrix)
  return (study.dat)
}

# read csv, convert dataframe to list of matrix
get_geno <- function(csv_path, sep=",") {
  
  # read csv-file with header into dataframe
  study.csv = read.csv2(file = csv_path , header = TRUE, sep=sep, dec=".")
  
  # get row names = first column of csv
  
  study.csv.data = study.csv[1:91,] # like cov.dat, traits.dat
  
  study.csv.rownames.list = as.list(study.csv.data[1])
  study.csv.rownames = study.csv.rownames.list[[1]]
  
  study.csv.data = study.csv.data[-c(1:6)] # remove first six columns FID-PHENOTYPE can be removed
  
  # convert strings to factors
  # study.sex_factor = as.factor(study.csv.data[2])
  # study.sex = as.integer(study.sex_factor)
  
  # not necessary for traits, no factors
  # sv = as.vector(study.csv.data[,2])
  # svi = as.integer(as.factor(sv))
  # study.csv.data[,2] = svi
  
  # convert to matrix
  nr = nrow(study.csv.data)
  nc = ncol(study.csv.data)
  study.matrix = as.matrix(study.csv.data, nrow=nr, ncol=nc)
  
  # remove all columns with NA values
  study.matrix.nona <- study.matrix[, !colSums(is.na(study.matrix))]
  
  # set rownames
  rownames(study.matrix.nona) = study.csv.rownames
  
  # return list of matrix
  study.dat = list(Study1=study.matrix.nona)
  return (study.dat)
}

# input file
# cov.f = "/home/kcan/UOC/tfm/mvgwas2-nf/data/mtar/cov.csv"
# traits.f = "/home/kcan/UOC/tfm/mvgwas2-nf/data/mtar/traits.csv"
# geno.f = "/home/kcan/UOC/tfm/mvgwas2-nf/data/mtar/geno.csv"

# cov.f = "/home/kcan/UOC/tfm/mvgwas2-nf/data/manta/eg.covariates.txt"
# cov_preproc.tsv
cov.f = "/home/kcan/UOC/tfm/mvgwas2-nf/data/mtar/cov_preproc.tsv"
traits.f = "/home/kcan/UOC/tfm/mvgwas2-nf/data/mtar/pheno_preproc.tsv"
geno.f ="/home/kcan/UOC/tfm/mvgwas2-nf/data/mtar/genotypes_plink1_rA.raw"

# convert input data to list of matrix
cov.dat = get_covariates(cov.f, sep="\t")
traits.dat = get_traits(traits.f, sep="\t")
geno.dat = get_geno(geno.f, sep=" ")

# traits:
#
# a numeric list, each sublist containing trait information for each study. In each study, a numeric 
# n×K matrix with each row as an independent individual (n) and each column as a separate trait (K). 
# If subject i is without trait k, the corresponding value is set as NA. The number of traits 
# in each study can be different but the names of traits are required.
#
# traits.dat:
##              LDL HDL TG
## SUBJ1 -2.6166166  NA NA
## SUBJ2  0.8519045  NA NA
## SUBJ3  1.2438145  NA NA

# covariates:
#
# a numeric list, each sublist containing covariates information for each study. 
# In each study, a numeric n×D matrix with each row as an independent individual (n) and 
# each column as a covariate (D).
#
# cov.dat:
##       cov1        cov2
## SUBJ1    1 -0.29713414
## SUBJ2    0 -0.08074092
## SUBJ3    0  0.55975760

# genotype:
#
# a numeric list, each sublist containing genotype information for each study. 
# In each study, a numeric n×m matrix with each row as an independent individual (n) and 
# each column as a SNP (m). Each genotype should be coded as 0, 1, 2, and NA for AA, Aa, aa, and missing, 
# where A and a represents a major and minor allele, respectively. 
# The number of SNPs in each study can be different but the names of SNPs are required. 
# Also, the number of studies must be the same in genotype, covariates and traits lists. 
# The order of subject ID must be the same among traits, covariates, and genotype within each study.
#
# geno.dat
##       SNP2148 SNP2149 SNP2150 SNP2151 SNP2152 SNP2153 SNP2154 SNP2155
## SUBJ1       0       0       0       0       0       0       0       0
## SUBJ2       0       0       0       0       0       0       0       0
## SUBJ3       0       0       0       0       0       0       0       0

# n number of individuals
# K number of traits
# D number of covariates
# m number of SNPs

# MTAR: Compute the summary statistics given the individual-level data
# returns: A list containing summary statistics for each trait. If covariance is TRUE, 
# the score summary statistics U and its covariance matrix V are returned.

obs.stat <- MTAR::Get_UV_from_data(traits = traits.dat,
                             covariates = cov.dat,
                             genotype = geno.dat,
                             covariance = TRUE)

# show results
# obs.stat$U$P1[1:3]
# obs.stat$V$P1[1:3]

U = obs.stat$U
V = obs.stat$V

# create MAF data (named vector)
maf_entries = length(U[[1]])
MAF = c(1:maf_entries)
MAF[1:maf_entries] = 0.012678 # just a sample value
# names(MAF) = c(paste0("SNP", 1:maf_entries))
names(MAF) = names(U[[1]])

# call MTAR
pval <-  MTAR(U = U, V = V, MAF = MAF)

# process results
print(pval)

# write pval values to a file
# pval$MTAR0, 
# pval$cMTAR$p, pval$cMTAR$rho1.min, pval$cMTAR$rho2.min
# pval$iMTAR$p, pval$iMTAR$rho1.min, pval$iMTAR$rho2.min
# pval$cctP



