// mvgwas2.nf

/*
 * Copyright (c) 2021, Diego Garrido-Martín
 *
 * This file is part of 'mvgwas-nf':
 * A Nextflow pipeline for multivariate GWAS using MANTA
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
 
 /* 2023, Kaan Candar: 
  * - upgade to Nextflow DSL2
  * - single Nextflow script for processing various multivariate GWAS methods:
  *   - MANTA
  *   - GEMMA
 */

// general parameters
debug_flag = false // set to true for additional logging

// workflow_conditional2.nf
// run with process provided on the command line:
// nextflow run workflow_conditional2.nf --process p1,p3

// proccess to be run
params.methods = null

// convert process string to a list
params.methodsList = params.methods?.split(',') as List

if (debug_flag)
{
    println "methodsList: " + params.methodsList
}

all_methods_list = ["manta","gemma","mtar","mostest"]


// Define common method parameters
params.help = false

params.geno = null
params.l = 500

params.pheno = null


// Define method specific parameters
// MANTA
params.cov = null
params.t = 'none'
params.i = 'none'
params.ng = 10
params.dir = 'result'
params.out = 'mvgwas.tsv'


// Print usage

if (params.help) {
    log.info ''
    log.info 'mvgwas2-nf: A pipeline for multivariate Genome-Wide Association Studies'
    log.info '========================================================================'
    log.info 'Performs multi-trait GWAS using various methods:'
    log.info ' - MANTA (https://github.com/dgarrimar/manta)'
    log.info ' - GEMMA'
    log.info ' - MTAR'
    log.info ' - MOSTEST'
    log.info ''
    log.info 'requires Nextflow DSL2'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run mvgwas2.nf [options]'
    log.info ''
    log.info 'Common parameters:'
    log.info ' --methods "manta,gemma,mtar,mostest" select one or more GWAS methods'
    log.info ' --pheno PHENOTYPES          phenotype file'
    log.info ' --geno GENOTYPES            indexed genotype VCF file'
    log.info " --l VARIANTS/CHUNK          variants tested per chunk (default: $params.l)"
    log.info ''
    log.info 'Parameters for MANTA:'
    log.info ' --cov COVARIATES            covariate file'
    log.info " --t TRANSFOMATION           phenotype transformation: none, sqrt, log (default: $params.t)"
    log.info " --i INTERACTION             test for interaction with a covariate: none, <covariate> (default: $params.i)"
    log.info " --ng INDIVIDUALS/GENOTYPE   minimum number of individuals per genotype group (default: $params.ng)"
    log.info " --dir DIRECTORY             output directory (default: $params.dir)"
    log.info " --out OUTPUT                output file (default: $params.out)"
    log.info ''
    exit(1)
}

// Check mandatory parameters
// parameters required for all methods:
if (params.methodsList == null)
{
    params.help
    exit 1, "GWAS method not specified."
}

if (! all_methods_list.containsAll(params.methodsList))
{
    params.help
    wrong_methods = params.methodsList.minus(all_methods_list)
    exit 1, "wrong GWAS methods specified: " + wrong_methods
}

if (!params.pheno) {
    params.help
    exit 1, "Phenotype file not specified."
} else if (!params.geno) {
    params.help
    exit 1, "Genotype not specified."
}

// Print parameter selection
log.info ''
log.info 'Parameters'
log.info '------------------'
log.info "Phenotype data               : ${params.pheno}"
log.info "Genotype data                : ${params.geno}"
log.info "Covariates                   : ${params.cov}"
log.info "Variants/chunk               : ${params.l}"
log.info "Phenotype transformation     : ${params.t}"
log.info "Interaction                  : ${params.i}"
log.info "Individuals/genotype         : ${params.ng}" 
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
log.info ''

// parameters required for specific methods
if ("manta" in params.methodsList) 
{
    if (!params.cov) {
        params.help
        exit 1, "Covariate file not specified for MANTA."
    }
}

// main workflow
workflow {

    println "methods to be processed: " + params.methods

    // input files
    fileGenoVcf = Channel.fromPath(params.geno)
    fileGenoTbi = Channel.fromPath("${params.geno}.tbi")
    filePheno = Channel.fromPath(params.pheno)

    fileCov = Channel.fromPath(params.cov)

    // common processing step
    common_p0()

    if ("split" in params.methodsList) {
        chunks = common_p0_split_genotype(fileGenoVcf, fileGenoTbi) | flatten
    }


    // specific processing steps for Manta
    if ("manta" in params.methodsList) {
        manta_p1()
        manta_p2()

        // preprocess phenotype and covariate data
        tuple_files = manta_p1_preprocess_pheno_cov(filePheno, fileCov, fileGenoVcf)
    }
    
    // specific processing steps for GEMMA
    if ("gemma" in params.methodsList) {
        gemma_p1()
        gemma_p2()
    }

    // specific processing steps for MTAR
    if ("mtar" in params.methodsList) {
        mtar_p1()
        mtar_p2()
    }

    // specific processing steps for MOSTEST
    if ("mostest" in params.methodsList) {
        mostest_p1()
        mostest_p2()
    }
}

process common_p0 {
    exec:
    println "this is process common_p0"
}

process manta_p1 {
    exec:
    println "this is process manta_p1"
}

process manta_p2 {
    exec:
    println "this is process manta_p2"
}

// Manta Step 1: Preprocess phenotype and covariate data
process manta_p1_preprocess_pheno_cov {
  
    debug debug_flag
    
    input:
    path pheno_file
    path cov_file
    path vcf_file
    
    output:
    tuple file("pheno_preproc.tsv.gz"), file("cov_preproc.tsv.gz")
    
    script:
    """
    echo "preprocessing files" $pheno_file $cov_file $vcf_file
    echo "preprocess.R --phenotypes $pheno_file --covariates $cov_file --genotypes $vcf_file --out_pheno pheno_preproc.tsv.gz --out_cov cov_preproc.tsv.gz --verbose"
    ${moduleDir}/bin/manta/preprocess.R --phenotypes $pheno_file --covariates $cov_file --genotypes $vcf_file --out_pheno pheno_preproc.tsv.gz --out_cov cov_preproc.tsv.gz --verbose
    """
}

process gemma_p1 {
    exec:
    println "this is process gemma_p1"
}

process gemma_p2 {
    exec:
    println "this is process gemma_p2"
}

process mtar_p1 {
    exec:
    println "this is process mtar_p1"
}

process mtar_p2 {
    exec:
    println "this is process mtar_p2"
}

process mostest_p1 {
    exec:
    println "this is process mostest_p1"
}

process mostest_p2 {
    exec:
    println "this is process mostest_p2"
}


// Split genotype VCF file
process common_p0_split_genotype {
  
  debug debug_flag 
  
  input:
  file(vcf) // from file(params.geno)
  file(index) // from file("${params.geno}.tbi")
  
  output:
  path("chunk*")
  
  script:

  if (debug_flag) {
    log.info "vcf-file: ${vcf}"  // e.g. eg.genotypes.vcf.gz
    log.info "tbi-file: ${index}"
    log.info "params.l: ${params.l}" // e.g. 500
  }
  
  """
  echo "splitting file" $vcf $index
  bcftools query -f '%CHROM\t%POS\n' $vcf > positions
  split -d -a 10 -l ${params.l} positions chunk
  """
}
