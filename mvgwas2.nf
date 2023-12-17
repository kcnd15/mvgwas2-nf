// mvgwas2.nf

/*
 * Copyright (c) 2021, Diego Garrido-Martín
 *
 * This file is part of 'mvgwas-nf':
 * A Nextflow pipeline for multivariate GWAS using various methods
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
params.debug_flag = false // set to true for additional logging

// workflow_conditional2.nf
// run with process provided on the command line:
// nextflow run workflow_conditional2.nf --process p1,p3

// proccess to be run
params.methods = null

// convert process string to a list
params.methodsList = params.methods?.split(',') as List

if (params.debug_flag)
{
    println "methodsList: " + params.methodsList
}

all_methods_list = ["manta","gemma","mtar","mostest"]


// Define common method parameters
params.help = false
params.l = 500 // No. of variants/chunk 

params.geno = null
params.pheno = null

params.dir = 'result' // output directory


// Define method specific parameters
// MANTA

params.cov = null
params.t = 'none' // phenotype transformation: none, sqrt, log
params.i = 'none' // test for interaction with a covariate: none, <covariate>
params.ng = 10 // minimum number of individuals per genotype group

params.manta_out = 'mvgwas.tsv'

// GEMMA
params.maf = 0.01 // MAF filter
params.gemma_out = 'gemma.tsv'
params.threads = 1 // No. of threads


// Print usage

if (params.help) {
    log.info ''
    log.info 'mvgwas2-nf: A pipeline for multivariate Genome-Wide Association Studies'
    log.info '========================================================================'
    log.info 'Performs multi-trait GWAS using various methods:'
    log.info ''
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
    log.info ' --pheno PHENOTYPES          (covariate-adjusted) phenotype file'
    log.info ' --geno GENOTYPES            indexed genotype VCF file'
    log.info " --dir DIRECTORY             output directory (default: $params.dir)"
    log.info " --l VARIANTS/CHUNK          variants tested per chunk (default: $params.l)"
    log.info ''
    log.info 'Parameters for MANTA:'
    log.info ' --cov COVARIATES            covariate file'
    log.info " --t TRANSFOMATION           phenotype transformation: none, sqrt, log (default: $params.t)"
    log.info " --i INTERACTION             test for interaction with a covariate: none, <covariate> (default: $params.i)"
    log.info " --ng INDIVIDUALS/GENOTYPE   minimum number of individuals per genotype group (default: $params.ng)"
    log.info " --manta_out OUTPUT          output file (default: $params.manta_out)"
    log.info ''
    log.info 'Parameters for GEMMA:'
    log.info " --maf MAF                   MAF filter (default: ${params.maf}"
    log.info " --threads THREADS           number of threads (default: ${params.threads})"
    log.info " --gemma_out OUTPUT          output file prefix (default: ${params.gemma_out})"
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

// parameters required for specific methods
if ("manta" in params.methodsList) 
{
    if (!params.cov) {
        params.help
        exit 1, "Covariate file not specified for MANTA."
    }
}

// Print parameter selection
log.info ''
log.info 'General parameters'
log.info '------------------'
log.info "GWAS methods                 : ${params.methodsList}"
log.info "Phenotype data               : ${params.pheno}"
log.info "Genotype data                : ${params.geno}"
log.info "Output directory             : ${params.dir}"
log.info ''

if ("manta" in params.methodsList) {

    log.info 'MANTA parameters'
    log.info '------------------'
    log.info "Covariates                   : ${params.cov}"
    log.info "Variants/chunk               : ${params.l}"
    log.info "Phenotype transformation     : ${params.t}"
    log.info "Interaction                  : ${params.i}"
    log.info "Individuals/genotype         : ${params.ng}" 
    log.info "Output file                  : ${params.manta_out}"
    log.info ''
}

if ("gemma" in params.methodsList) {

    log.info 'GEMMA parameters'
    log.info '------------------'
    log.info "MAF                          : ${params.maf}"
    log.info "Threads                      : ${params.threads}"
    log.info "Output file                  : ${params.gemma_out}"
    log.info ''
}

// main workflow
workflow {

    // input files
    fileGenoVcf = Channel.fromPath(params.geno)
    fileGenoTbi = Channel.fromPath("${params.geno}.tbi")
    filePheno = Channel.fromPath(params.pheno)

    fileCov = Channel.fromPath(params.cov)

    // common processing step
    chunks = common_p0_split_genotype(fileGenoVcf, fileGenoTbi) | flatten

    // specific processing steps for Manta
    if ("manta" in params.methodsList) {

        // preprocess phenotype and covariate data
        tuple_files = manta_p1_preprocess_pheno_cov(filePheno, fileCov, fileGenoVcf)

        // perform multivariate GWAS analysis using MANTA and output results
        manta_p2_mvgwas_manta(tuple_files, fileGenoVcf, fileGenoTbi, chunks) | manta_p3_collect_summary_statistics
    }
    
    // specific processing steps for GEMMA
    if ("gemma" in params.methodsList) {

        // Preprocess genotype and phenotype data
        tuple_files = gemma_p1_preprocess_geno_pheno(fileGenoVcf, filePheno)

        // Compute kinship
        tuple_eigen = gemma_p2_kinship(tuple_files)

        // // Test for association between phenotypes and genetic variants
        // and collect summary statistics
        gemma_p3_test_association(tuple_files, tuple_eigen, chunks) | gemma_p4_collect_summary_statistics
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

// ------------------------------------------------------------------------------
// common processing
// ------------------------------------------------------------------------------

// Split genotype VCF file
process common_p0_split_genotype {
  
  debug params.debug_flag 
  
  input:
  file(vcf) // from file(params.geno)
  file(index) // from file("${params.geno}.tbi")
  
  output:
  path("chunk*")
  
  script:

  if (params.debug_flag) {
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

// ------------------------------------------------------------------------------
// MANTA processing
// ------------------------------------------------------------------------------

// MANTA Step 1: Preprocess phenotype and covariate data
process manta_p1_preprocess_pheno_cov {
  
    debug params.debug_flag
    
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

// MANTA Step 2: Test for association between phenotypes and genetic variants using MANTA
process manta_p2_mvgwas_manta {
  
    debug params.debug_flag 

    input:

    tuple file(pheno), file(cov) // from preproc_ch
    file(vcf) // from file(params.geno)
    file(index) // from file("${params.geno}.tbi")
    each path(chunk) // file(chunk) // from chunks_ch

    output:

    path('sstats.*.txt') // optional true // into sstats_ch
    
    script:
    
    if (params.debug_flag) {
      log.info "logging processing of file ${chunk}"
    }
    
    """
    echo "processing file " $chunk
    
    chunknb=\$(basename $chunk | sed 's/chunk//')
    echo "chunknb:" \$chunknb
    
    # check for number of chromosomes
    if [[ \$(cut -f1 $chunk | sort | uniq -c | wc -l) -ge 2 ]]; then
        echo "entering if..."
        k=1
        cut -f1 $chunk | sort | uniq | while read chr; do
          region=\$(paste <(grep -P "^\$chr\t" $chunk | head -n 1) <(grep -P "^\$chr\t" $chunk | tail -n 1 | cut -f2) | sed 's/\t/:/' | sed 's/\t/-/')
        
          echo "if region: [" \$region "]"
          echo "--output" sstats.\$k.tmp
          echo "test.R --phenotypes $pheno --covariates $cov --genotypes $vcf --region "\$region" --output sstats.\$k.tmp --min_nb_ind_geno ${params.ng} -t ${params.t} -i ${params.i} --verbose"
        
          ${moduleDir}/bin/manta/test.R --phenotypes $pheno --covariates $cov --genotypes $vcf --region "\$region" --output sstats.\$k.tmp --min_nb_ind_geno ${params.ng} -t ${params.t} -i ${params.i} --verbose 
        
          ((k++))
      done
      cat sstats.*.tmp > sstats.\${chunknb}.txt
    else
        # only 1 chromosome
        echo "entering else..."
        
        region=\$(paste <(head -n 1 $chunk) <(tail -n 1 $chunk | cut -f2) | sed 's/\t/:/' | sed 's/\t/-/')
        
        echo "else region: [" \$region "]"
        echo "--output" sstats.\${chunknb}.txt
        echo "test.R --phenotypes $pheno --covariates $cov --genotypes $vcf --region "\$region" --output sstats.\${chunknb}.txt --min_nb_ind_geno ${params.ng} -t ${params.t} -i ${params.i} --verbose"
        
        ${moduleDir}/bin/manta/test.R --phenotypes $pheno --covariates $cov --genotypes $vcf --region "\$region" --output sstats.\${chunknb}.txt --min_nb_ind_geno ${params.ng} -t ${params.t} -i ${params.i} --verbose
    fi
    """
}

// MANTA step 3: collect resulting summary statistics
process manta_p3_collect_summary_statistics {
  
    // creates an output text file containing the multi-trait GWAS summary statistics

    debug params.debug_flag
    publishDir "${params.dir}", mode: 'copy'     

    input:
    file(out) // from pub_ch

    output:
    file(out) // into end_ch

    script:
    
    if (params.debug_flag) {
      log.info("process end")
      log.info("params.i: ${params.i}")
      log.info("input file: ${out}")
    }
    
    if (params.i == 'none')
    """
    sed -i "1 s/^/CHR\tPOS\tID\tREF\tALT\tF\tR2\tP\\n/" ${out}
    """
    else
    """
    sed -i "1 s/^/CHR\tPOS\tID\tREF\tALT\tF($params.i)\tF(GT)\tF(${params.i}:GT)\tR2($params.i)\tR2(GT)\tR2(${params.i}:GT)\tP($params.i)\tP(GT)\tP(${params.i}:GT)\\n/" ${out}
    """
}


// ------------------------------------------------------------------------------
// GEMMA processing
// ------------------------------------------------------------------------------

// GEMMA Step 1: Preprocess genotype and phenotype data
process gemma_p1_preprocess_geno_pheno {

    cpus params.threads

    input:
    path vcf // file vcf from file(params.geno)    
    path pheno // file pheno from file(params.pheno)

    output:
    tuple file("geno.bed"), file("geno.bim"), file("geno.fam") // into geno_ch

    script:
    """
    comm -12 <(bcftools view -h $vcf | grep CHROM | sed 's/\\t/\\n/g' | sed '1,9d' | sort) <(cut -f1 $pheno | sort) > keep.txt
    plink2 --vcf $vcf --keep keep.txt --make-bed --out geno --threads ${params.threads} 
    awk '{print int(NR)"\t"\$0}' <(cut -f2 geno.fam) > idx
    join -t \$'\t' -1 2 -2 1 <(sort -k2,2 idx) <(sort -k1,1 $pheno) | sort -k2,2n | cut -f1,2 --complement > pheno.tmp
    paste <(cut -f1-5 geno.fam) pheno.tmp > tmpfile; mv tmpfile geno.fam
    """
}

// GEMMA Step 2: Compute kinship, obtain kinship matrix
process gemma_p2_kinship {

    cpus params.threads

    input:
    tuple file(bed), file(bim), file(fam) // from geno_ch

    output:
    tuple file("kinship.sXX.eigenD.txt"), file("kinship.sXX.eigenU.txt") // into kinship_ch

    script:
    """
    # Compute kinship
    export OPENBLAS_NUM_THREADS=${params.threads}
    prefix=\$(basename $bed | sed 's/.bed//')
    plink2 --bfile \$prefix --maf ${params.maf} --indep-pairwise 50 5 0.8 --threads ${params.threads}
    plink2 --bfile \$prefix --extract plink2.prune.in --out geno.pruned --make-bed --threads ${params.threads}
    gemma -gk 2 -bfile geno.pruned -outdir . -o kinship
    gemma -bfile geno.pruned -k kinship.sXX.txt -eigen -outdir . -o kinship.sXX
    """
}

// GWAS: testing (GEMMA)
// GEMMA Step 3: Test for association between phenotypes and genetic variants

process gemma_p3_test_association {

    cpus params.threads

    input:
    tuple file(bed), file(bim), file(fam) // from geno_ch
    tuple file(kinship_d), file(kinship_u) // from kinship_ch
    each file(chunk) // from chunks_ch

    output:
    file('gemma.0*.assoc.txt') // into sstats_ch

    script:
    """
    export OPENBLAS_NUM_THREADS=${params.threads}
    pids=\$(seq \$(echo \$(awk '{print NF}' $fam | head -1) - 5 | bc -l))
    chunknb=\$(basename $chunk | sed 's/chunk//')

    if [[ \$(cut -f1 $chunk | sort | uniq -c | wc -l) -ge 2 ]]; then
        k=1
        cut -f1 $chunk | sort | uniq | while read chr; do
            paste <(grep -P "^\$chr\t" $chunk | head -1) <(grep -P "^\$chr\t" $chunk | tail -1 | cut -f2) > region
            plink2 -bfile geno --extract bed1 region --make-bed --out geno.ss --threads ${params.threads}
            paste <(cut -f1-5 geno.ss.fam) <(cut -f1-5 --complement geno.fam) > tmpfile; mv tmpfile geno.ss.fam
            (timeout 120 gemma -lmm -b geno.ss -d $kinship_d -u $kinship_u -n \$pids -outdir . -o gemma.k\$k -maf ${params.maf} &> STATUS || exit 0)
            if [[ \$(grep ERROR STATUS) ]]; then
                touch gemma.k\$k.assoc.txt
            else
                gemma -lmm -b geno.ss -d $kinship_d -u $kinship_u -n \$pids -outdir . -o gemma.k\$k -maf ${params.maf}
            fi
            ((k++))
        done
        cat gemma.k*.assoc.txt > gemma.\${chunknb}.assoc.txt        
    else
        paste <(head -1 $chunk) <(tail -1 $chunk | cut -f2) > region
        plink2 -bfile geno --extract bed1 region --make-bed --out geno.ss --threads ${params.threads}
        paste <(cut -f1-5 geno.ss.fam) <(cut -f1-5 --complement geno.fam) > tmpfile; mv tmpfile geno.ss.fam
        (timeout 120 gemma -lmm -b geno.ss -d $kinship_d -u $kinship_u -n \$pids -outdir . -o gemma.\${chunknb} -maf ${params.maf} &> STATUS || exit 0)
        if [[ \$(grep ERROR STATUS) ]]; then
            touch gemma.\${chunknb}.assoc.txt
        else
            gemma -lmm -b geno.ss -d $kinship_d -u $kinship_u -n \$pids -outdir . -o gemma.\${chunknb} -maf ${params.maf}
        fi
    fi
    """
}

// GEMMA Step 4: Collect summary statistics

process gemma_p4_collect_summary_statistics {

    publishDir "${params.dir}", mode: 'copy'

    input:
    file(out) // from pub_ch

    output:
    file(out) // into end_ch

    script:
    """
    head -1 $out > header
    grep -v n_miss $out > tmpfile; cat header tmpfile > $out
    """
}

// ------------------------------------------------------------------------------
// MTAR processing
// ------------------------------------------------------------------------------

process mtar_p1 {
    exec:
    println "this is process mtar_p1"
}

process mtar_p2 {
    exec:
    println "this is process mtar_p2"
}

// ------------------------------------------------------------------------------
// MOSTest processing
// ------------------------------------------------------------------------------

process mostest_p1 {
    exec:
    println "this is process mostest_p1"
}

process mostest_p2 {
    exec:
    println "this is process mostest_p2"
}

