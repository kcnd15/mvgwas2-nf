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


// workflow_conditional2.nf
// run with process provided on the command line:
// nextflow run workflow_conditional2.nf --process p1,p3

// proccess to be run
params.methods = "manta,gemma"

// convert process string to a list
params.methodsList = params.methods?.split(',') as List

println "this is main"
println "methodsList: " + params.methodsList

// general parameters
debug_flag = false // set to true for additional logging

// Define common method parameters
params.geno = null
params.l = 500

// Define method specific parameters
// ...


// main workflow
workflow {

    println "this is the main workflow"
    println "methods: " + params.methods

    // input files
    fileGenoVcf = Channel.fromPath(params.geno)
    fileGenoTbi = Channel.fromPath("${params.geno}.tbi")

    // common processing step
    common_p0()

    if ("split" in params.methodsList) {
        chunks = common_p0_split_genotype(fileGenoVcf, fileGenoTbi) | flatten
    }


    // specific processing steps for Manta
    if ("manta" in params.methodsList) {
        manta_p1()
        manta_p2()
    }
    
    if ("gemma" in params.methodsList) {
        gemma_p1()
        gemma_p2()
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

process gemma_p1 {
    exec:
    println "this is process gemma_p1"
}

process gemma_p2 {
    exec:
    println "this is process gemma_p2"
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
