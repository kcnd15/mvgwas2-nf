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
  *   - MTAR
  *   - MOSTest
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
params.result_dir = null


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

// MTAR
params.mtar_cov = null
params.datadir = '.'

// MOSTest
params.mostest_pheno = null
params.mostest_bfile = null
params.mostest_out_prefix = null
params.mostest_data_dir = null
params.mostest_result_dir = null
params.mostest_sample = null


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
    log.info " --maf MAF                   MAF filter (default: ${params.maf})"
    log.info " --threads THREADS           number of threads (default: ${params.threads})"
    log.info " --gemma_out OUTPUT          output file prefix (default: ${params.gemma_out})"
    log.info ''
    log.info 'Parameters for MTAR:'
    log.info " --cov COVARIATES            covariate file"
    log.info " --datadir dir               data directory (default: ${params.datadir})"
    log.info ''
    log.info 'Parameters for MOSTest:'
    log.info " --mostest_pheno pheno       phenotype file (default: ${params.mostest_pheno})"
    log.info " --mostest_bfile bfile       PLINK bfile prefix (default: ${params.mostest_bfile})"
    log.info " --mostest_out_prefix p      result files prefix (default: ${params.mostest_out_prefix})"
    log.info " --mostest_data_dir dir      data directory (default: ${params.mostest_data_dir})"
    log.info " --mostest_result_dir dir    result directory (default: ${params.mostest_result_dir})"
    log.info " --mostest_sample            use MOSTest sample data (default: ${params.mostest_sample})"
  
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

if ("mostest" in params.methodsList) {

    log.info 'MOSTest parameters'
    log.info '--------------------'
    log.info "Phenotype file               : ${params.mostest_pheno}"
    log.info "PLINK bfile prefix           : ${params.mostest_bfile}"
    log.info "Output file prefix           : ${params.mostest_out_prefix}"
    log.info "Input data directory         : ${params.mostest_data_dir}"
    log.info "Result directory             : ${params.mostest_result_dir}"
    log.info ''
}

// main workflow
workflow {

    // input files
    fileGenoVcf = Channel.fromPath(params.geno)
    fileGenoTbi = Channel.fromPath("${params.geno}.tbi")
    filePheno = Channel.fromPath(params.pheno)

    fileCov = Channel.fromPath(params.cov)
    
    fileMostestPheno = Channel.fromPath(params.mostest_pheno)

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

        // Test for association between phenotypes and genetic variants
        // and collect summary statistics
        gemma_p3_test_association(tuple_files, tuple_eigen, chunks) | gemma_p4_collect_summary_statistics
    }

    // specific processing steps for MTAR
    if ("mtar" in params.methodsList) {
      
        // preprocess phenotype and covariate data
        tuple_files = manta_p1_preprocess_pheno_cov(filePheno, fileCov, fileGenoVcf)
      
        mtar_p1_prepare_summary_statistics(filePheno, fileCov, fileGenoVcf)

        mtar_p2_calculate_cov()
        mtar_p3_run_mtar()
    }

    // specific processing steps for MOSTEST
    if ("mostest" in params.methodsList) {
        if (params.geno) {
          (bfile, tuple_bfiles) = mostest_create_plink_files(fileGenoVcf)
          log.info("tuple: $bfile, $tuple_bfiles")
        } else {
          bfile = params.mostest_bfile
        }
        
        out_prefix = mostest_run_mostest(fileMostestPheno, bfile, tuple_bfiles)
        mostest_process_results(bfile, out_prefix)
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

// MTAR Step 1: Prepare Summary Statistics for Multiple Traits
process mtar_p1_prepare_summary_statistics {
  
    debug params.debug_flag
    
    input:
    path pheno_file
    path cov_file
    path vcf_file

    script:
    
    log.info("mtar_p1_prepare_summary_statistics - MTAR Step 1: Prepare Summary Statistics for Multiple Traits")
    log.info("pheno_file: ${pheno_file}")
    log.info("cov_file: ${cov_file}")
    log.info("vcf_file: ${vcf_file}")
    
    """
    # generate PLINK 1.9 files
    plink --vcf ${vcf_file} --out genotypes_plink1
    
    # generate Genotype matrix
    plink --bfile genotypes_plink1 --recode A --out genotypes_plink1_rA

    # MTAR prepare summary statistics
    ${moduleDir}/bin/mtar/prepareSumStat.R -d ${params.datadir} -c $cov_file -t $pheno_file -g genotypes_plink1_rA.raw
    """
}

// MTAR Step 2: Calculate Covariances of Z-scores between Traits from Overlapping Samples
process mtar_p2_calculate_cov {
    exec:
    println "Step 2: Calculate Covariances of Z-scores"
}

// MTAR Step 3. Run Multi-Traits Rare-Variant Association Tests (MTAR)
process mtar_p3_run_mtar {
    exec:
    println "Step 3. Run MTAR"
}

// ------------------------------------------------------------------------------
// MOSTest processing
// ------------------------------------------------------------------------------

// params.mostest_pheno = 'pheno.txt' 
// params.mostest_bfile = 'chr21'
// params.mostest_out_prefix = 'mostest_results'
// params.mostest_data_dir = '$baseDir/data/mostest'
// params.mostest_result_dir = 'result'

process mostest_create_plink_files {
  
    debug params.debug_flag
    
    input:
    path vcf_file
    
    output:
    val "genotypes_plink1"
    tuple file("genotypes_plink1.bed"), 
          file("genotypes_plink1.bim"),
          file("genotypes_plink1.fam")
    
    script:
    println "mostest_create_plink_files: vcf_file = $vcf_file"
    
    """
    # generate PLINK 1.9 files
    echo "mostest_create_plink_files"
    pwd
    plink --vcf ${vcf_file} --out genotypes_plink1
    """
}

process mostest_run_mostest {
  
    debug params.debug_flag
    
    input:
      file(pheno)
      val(bfile)
      tuple file("genotypes_plink1.bed"), 
          file("genotypes_plink1.bim"),
          file("genotypes_plink1.fam")
      
    output:
      // tuple file("${params.mostest_out_prefix}.mat"), file("${params.mostest_out_prefix}_zmat.mat")
      val "${params.mostest_out_prefix}"
    
    script:
    println "this is process mostest_run_mostest: bfile = $bfile"
    println "prefix = ${params.mostest_out_prefix}"
    
    """
    echo "mostest_run_mostest: pheno= $pheno"
    pwd
    # cd ${moduleDir}/bin/mostest
    export MATLABPATH=/home/j/Documents/MATLAB/mine
    # echo "MOSTest pheno: ${params.mostest_data_dir}"
    
    export MATLABPATH=${moduleDir}
    echo "mostest_run_mostest, MATLABPATH:"
    
    # ${moduleDir}/bin/mostest/run_mostest.sh $pheno ${params.mostest_bfile} ${params.mostest_out_prefix} ${params.mostest_data_dir} ${params.mostest_result_dir}
    ${moduleDir}/bin/mostest/run_mostest.sh $pheno $bfile ${params.mostest_out_prefix} . ${params.mostest_result_dir}
    
    # ln -s ${params.mostest_result_dir}/${params.mostest_out_prefix}.mat .
    # ln -s ${params.mostest_result_dir}/${params.mostest_out_prefix}_zmat.mat .
    """
}

process mostest_process_results {
  
    debug params.debug_flag
    
    input:
    // tuple file(result_mat), file(result_zmat)
    val bfile
    val out_prefix
  
    script:
    println "this is process mostest_process_results:"
    println "bfile: $bfile, out_prefix: $out_prefix"
    
    """
    echo ${moduleDir}/bin/mostest/process_results.py ${params.data_dir}/${bfile}.bim ${params.result_dir}/$out_prefix
        
    python3 ${moduleDir}/bin/mostest/process_results.py ${params.data_dir}/${bfile}.bim ${params.result_dir}/$out_prefix

    # python3 bin/mostest/process_results.py data/mostest/chr21.bim result/mostest_results -> produces an error
    """
}

