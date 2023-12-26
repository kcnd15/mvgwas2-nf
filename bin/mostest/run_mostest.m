% To run MOSTest, set your current folder to <MOST_ROOT>. 
% Then start matlab, define pheno, out, bfile, snps and nsubj variables as 
% shown below, and execute mostest.m script:


function exitcode = run_mostest(pheno_file, bfile_prefix, out_prefix, data_dir, result_dir)

  fprintf("starting mostest.m .....\n\n")

  pheno = pheno_file;            % full or relative path to the phenotype file; 'pheno.txt'
  bfile = bfile_prefix;          % full or relative path to plink bfile prefix; 'chr21'
  out = out_prefix;              % prefix for the output files; 'results'
  
  fprintf("phenotype file      : %s\n", pheno)
  fprintf("plink bfile prefix  : %s\n", bfile)
  fprintf("output files prefix : %s\n", out)
  fprintf("data directory      : %s\n", data_dir)
  fprintf("result directory    : %s\n", data_dir)
  fprintf("\n")
  
  mostest(pheno, bfile, out, data_dir, result_dir); % start the MOSTest analysis
        
  exitcode = 0;
end