function mostest(pheno_file, bfile_prefix, out_prefix, data_dir, result_dir)

% Phenotype file               : pheno_file, eg. mvgwas2-nf/data/eg.phenotypes.txt
% PLINK bfile prefix           : bfile_prefix, e.g. chr21
% Output file prefix           : out_prefix, e.g. mostest_results_chr21
% Input data directory         : data_dir, e.g. mvgwas2-nf/data/mostest
% Result directory             : result_dir, e.g. mvgwas2-nf/result

%  phenotype file      : phenotypes.clean.tsv
%  plink bfile prefix  : genotypes_plink1
%  output files prefix : mostest_results_genotypes
%  data directory      : .
%  result directory    : /home/kcan/UOC/tfm/mvgwas2-nf/result/ADNI
%  
%  mostest.m: if isempty(zmat_name)
%  trying to open file ./genotypes_plink1.bim
%  ...opened
%  trying to open file ./genotypes_plink1.fam
%  ...opened
%  500 snps and 808 subjects detected in bfile
%  Loading phenotype matrix from ./phenotypes.clean.tsv...
%  run_mostest.sh: MATLAB:table:RowIndexOutOfRange: Row index exceeds table dimensions.
%
% Work dir:
%  /home/kcan/UOC/tfm/mvgwas2-nf/work/92/b5551f1555c7e222fa5fbb5f3e2f6d

% direct settings for debugging one processing step in a 
% Nextflow working directory:
debug_workdir = false;
if debug_workdir
    pheno_file = "phenotypes.clean.tsv";
    bfile_prefix = "genotypes_plink1";
    out_prefix = "mostest_results_genotypes";
    data_dir = ".";
    result_dir = "/home/kcan/UOC/tfm/mvgwas2-nf/result/ADNI";
end

debug_flag = false;

if debug_flag
  fprintf("mostest.m: started...\n")
  if batchStartupOptionUsed
      fprintf("batchStartupOptionUsed\n")
  end
  
  % license check
  license_code = license('test','Statistics_Toolbox');
  fprintf("license('test','Statistics_Toolbox'): %d\n", license_code)

  fprintf("pheno_file   : %s\n", pheno_file)
  fprintf("bfile_prefix : %s\n", bfile_prefix)
  fprintf("out_prefix   : %s\n", out_prefix)
  fprintf("data_dir     : %s\n", data_dir)
  fprintf("result_dir   : %s\n", result_dir)
  fprintf("\n")
end

out = out_prefix;
pheno_file_name = pheno_file;
bfile = bfile_prefix;

% =============== parameters section =============== 

% optional arguments
if ~exist('zmat_name', 'var'), zmat_name = ''; end;                               % re-use univariate GWAS results from previous MOSTest analysis
if ~exist('chunk', 'var'), chunk = 10000; end;                                    % chunk size (how many SNPs to read at a time)
if ~exist('num_eigval_to_keep', 'var'), num_eigval_to_keep = 0; end;              % how many largest eigenvalues of C0 matrix (z score correlation) to keep, the remaining will be assigned to the num_eigval_to_keep-th eigenvalue, num_eigval_to_keep = 0 - keep all
if ~exist('apply_int', 'var'), apply_int = true; end;                             % apply rank-based inverse normal transform
if ~exist('use_pheno_corr', 'var'), use_pheno_corr = false; end;                  % use correlation structure of the phenotypes
if ~exist('auto_compile_shuffle', 'var'), auto_compile_shuffle = 1; end;          % automatically compile shuffle.mex
if ~exist('use_paretotails', 'var'), use_paretotails = false; end;                % use paretotails instead of the gamma and beta functions to fit the distribution of the MOSTest & minP test statistic under null
if ~exist('maf_threshold', 'var'), maf_threshold = 0.005; end;                      % ignore all variants with maf < maf_threshold in MOSTest analysis

% required input
if ~exist('out', 'var'),   error('out file prefix is required'); end
if isempty(zmat_name)
  if ~exist('pheno_file_name', 'var'), error('pheno file is required'); end
  if ~exist('bfile', 'var'), error('bfile is required'); end
end

% debug features - internal use only
if ~exist('perform_cca', 'var'), perform_cca = false; end;  % perform canonical correlation analysis
if ~exist('lam_reg', 'var'), lam_reg = nan; end;  %  default is to disable pre-whitening filter
if ~exist('snps', 'var'), snps=nan; end;                                          % number of SNPs in the analysis
if ~exist('nsubj', 'var'), nsubj=nan; end;                                        % number of subjects in the analysis
if ~exist('paretotails_quantile', 'var'), paretotails_quantile = 0.9999; end;       % a number close to 1.0, used as a second argument in MATLAB's paretotails
      
% =============== end of parameters section =============== 

if debug_flag
  fprintf("mostest.m: shuffle start...\n")
end

if auto_compile_shuffle && (exist('Shuffle') ~= 3), mex 'Shuffle.c'; end;   % ensure Shuffle is compiled

if debug_flag
  fprintf("mostest.m: shuffle end.\n")
end

tic

if isempty(zmat_name)

  fprintf("mostest.m: if isempty(zmat_name)\n")

  try
      bfile_bim = sprintf('%s/%s.bim', data_dir, bfile);
      fprintf("trying to open file %s\n", bfile_bim)
      fileID = fopen(bfile_bim);
      fprintf("...opened\n")
      bim_file = textscan(fileID,'%s %s %s %s %s %s');
      fclose(fileID);
  catch me
      % On error, print error message and exit with failure
      fprintf('mostest.m: current directory is %s\n', pwd)
      fprintf('mostest.m: ls %s\n', ls)
      fprintf('mostest.m %s / %s for file %s\n', me.identifier, me.message, bfile_bim)
      exit(1)
  end


  if isfinite(snps) && (snps ~= length(bim_file{1})), error('snps=%i is incompatible with .bim file; please check your snps parameter (or remove it to auto-detect #snps)', snps);end
  snps=length(bim_file{1});

  try
    bfile_fam = sprintf('%s/%s.fam', data_dir, bfile);
    fprintf("trying to open file %s\n", bfile_fam)
    fileID = fopen(bfile_fam);
    fprintf("...opened\n")
    fam_file = textscan(fileID,'%s %s %s %s %s %s');
    fclose(fileID);
  catch me
      % On error, print error message and exit with failure
      fprintf('mostest.m %s / %s for file %s\n', me.identifier, me.message, bfile_fam)
      exit(1)
  end
  
  if isfinite(nsubj) && (nsubj ~= length(fam_file{1})), error('nsubj=%i is incompatible with .fam file; please check your snps parameter (or remove it to auto-detect nsubj)', nsubj);end
  nsubj=length(fam_file{1});

  fprintf('%i snps and %i subjects detected in bfile\n', snps, nsubj);
  
  pheno = sprintf('%s/%s', data_dir, pheno_file_name);

  fprintf('Loading phenotype matrix from %s...\n', pheno);
  if 1 
      ymat_df = readtable(pheno, 'FileType', 'text', 'Delimiter', 'tab');
      % t = readtable("data.tsv", "FileType","text",'Delimiter', '\t');
      measures = ymat_df.Properties.VariableNames;
      
      % kc: remove first column "ID" if present 
      if measures{1} == 'ID'

        % ymat must have the same number of subjects as the .FAM file
        % ADNI data error: Row index exceeds table dimensions.
        keep_ids = fam_file{2}; % subject IDs in FAM-file, e.g. 385

        % take only those entries of the pheno-data in ymat_df,
        % whose IDs (column "ID", e.g. '002_S_0413') are also present in 
        % the fam-file (fam_file{1,2}) with ID column, e.g. '037_S_0501'

        % keep_ids 385x1 cell
        % keep_ids = {'002_S_0413', '002_S_1155', '002_S_4229'};
        % keep_ids_arr = table();
        % keep_ids_arr.ID = vertcat(keep_ids{:});
        % ymat_df(1,1) and keep_ids{1}, idx must be n*1 logical
        % idx = any(keep_ids_arr == ymat_df.ID(1));
        
        count_keep_ids = height(keep_ids);
        count_ymat_df = height(ymat_df);
        keep_idx = false(count_ymat_df,1);
        
        for row_ymat_df = 1:count_ymat_df
        
            for row_keep_ids = 1:count_keep_ids
                ymat_df_ID = ymat_df{row_ymat_df,1}{1};
                keep_ids_ID = keep_ids{row_keep_ids};
                if ymat_df_ID == keep_ids_ID
                    keep_idx(row_ymat_df) = true;
                    break
                end
            end
        end

        ymat_keep_df = ymat_df(keep_idx,:);

        % then remove the ID-column
        measures_length = length(measures);
        ymat_df =  ymat_keep_df(:,2:measures_length);

        % continue with normal processing
        % pheno_entries = height(ymat_df);
        % min_subj_pheno = min(pheno_entries, nsubj);
        % ymat_df = ymat_df(1:min_subj_pheno,:);
      end

      % continue with reduced table
      ymat_orig = table2array(ymat_df);
  else
      % an alternative helper code that reads phenotype matrix without a header
      ymat_orig=dlmread(pheno); ymat_orig=ymat_orig(:, 2:end);
      measures = cell(size(ymat_orig, 2), 1);
      for i=1:length(measures), measures{i} = sprintf('V%i', i); end;
  end
  npheno=size(ymat_orig, 2);
  fprintf('Done, %i phenotypes found\n', npheno);
  if size(ymat_orig, 1) ~= nsubj, error('roi matrix has info for %i subjects, while nsubj argument is specified as %i. These must be consistent.', size(ymat_orig, 1), nsubj); end;

  keep = (min(ymat_orig)~=max(ymat_orig));
  fprintf('Remove %i phenotypes (no variation)\n', length(keep) - sum(keep));
  ymat_orig = ymat_orig(:, keep);
  measures = measures(keep);
  npheno=size(ymat_orig, 2);

  % perform rank-based inverse normal transform, equivalently to the following R code:
  % DM[,f] <- qnorm(ecdf(DM[,f])(DM[,f]) - 0.5/dim(DM)[1])
  kurt = nan(npheno, 2);
  for pheno_index=1:npheno
    vals = ymat_orig(:, pheno_index); kurt(pheno_index, 1) = kurtosis(vals);
    if apply_int
      [F, X] = ecdf(vals); F=transpose(F(2:end)); X=transpose(X(2:end));
      vals = norminv(interp1(X,F,vals,'nearest') - 0.5 / length(vals));
    end
    ymat_orig(:, pheno_index) = vals; kurt(pheno_index, 2) = kurtosis(vals);
  end
  fprintf('kurtosis before INT - %.2f %.2f (mean, max)\n', mean(kurt(:, 1)), max(kurt(:, 1)))
  if apply_int, fprintf('kurtosis after  INT - %.2f %.2f (mean, max)\n', mean(kurt(:, 2)), max(kurt(:, 2))); end;

  if isfinite(lam_reg)
    C = corr(ymat_orig);
    C_reg = (1-lam_reg)*C + lam_reg*diag(max(0.01,diag(C))); % Ridge regularized covariance matrix
    C_inv = inv(C_reg);
    W_wht = chol(C_inv); % Whitening filter
    ymat = ymat_orig*W_wht'; % Whitened residualized data
  else
    ymat = ymat_orig;
  end

  if use_pheno_corr
    ymat_corr = corr(ymat);
  else
    ymat_corr = 'not computed';
  end

  fprintf('Perform GWAS on %s (%i SNPs are expected)...\n', bfile, snps)
  zmat_orig=zeros(snps, npheno, 'single');
  zmat_perm=zeros(snps, npheno, 'single');
  beta_orig=zeros(snps, npheno, 'single');  % skip saving p-values and standard errors (SE)
  beta_perm=zeros(snps, npheno, 'single');  % (can be derived from Z and BETA)
  nvec=zeros(snps, 1, 'single');
  freqvec=zeros(snps, 1, 'single');
  zvec_cca=nan(snps, 2);

  for i=1:chunk:snps
    j=min(i+chunk-1, snps);
    fprintf('gwas: loading snps %i to %i... ', i, j);    tic;
    
    try
        bfile_path = sprintf('%s/%s', data_dir, bfile);
        
        if debug_flag
          fprintf("\nPlinkRead_binary2, bfile=%s\n", bfile_path)
        end
        
        geno_int8 = PlinkRead_binary2(nsubj, i:j, bfile_path);
    catch me
        % On error, print error message and exit with failure
        fprintf('mostest.m %s / %s for file %s\n', me.identifier, me.message, bfile_path)
        exit(1)
    end
    
    fprintf('processing... ', i, j);   
    geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;

    shuffle_geno = Shuffle(geno);
    [rmat_orig_chunk, zmat_orig_chunk] = nancorr(ymat, geno);
    [rmat_perm_chunk, zmat_perm_chunk] = nancorr(ymat, shuffle_geno);

    if perform_cca
      fprintf('cca... ');   
      ymat1 = [ymat, ones(size(ymat, 1), 1)];
      for k=i:j
        % These two are equivalent:
        % [b, bint, r, rint, stats] = regress(y, [X ones(n, 1)]);  stats(3)      % based on F-test
        % [A, B, r, U, V, statsCCA] = canoncorr(X, y);             statsCCA.p  
        [b, bint, r, rint, stats] = regress(geno(:,         k-i+1), ymat1); zvec_cca(k, 1) = stats(3);
        [b, bint, r, rint, stats] = regress(shuffle_geno(:, k-i+1), ymat1); zvec_cca(k, 2) = stats(3);
      end
    end

    zmat_orig(i:j, :) = zmat_orig_chunk';
    zmat_perm(i:j, :) = zmat_perm_chunk';
    
    % https://stats.stackexchange.com/questions/32464/how-does-the-correlation-coefficient-differ-from-regression-slope
    beta_factor = std(ymat)' * (1./std(geno, 'omitnan'));
    beta_orig(i:j, :) = transpose(rmat_orig_chunk .* beta_factor);
    beta_perm(i:j, :) = transpose(rmat_perm_chunk .* beta_factor);
    nvec(i:j) = sum(isfinite(geno))';
    freqvec(i:j) = (1*sum(geno==1) + 2*sum(geno==2))' ./ (2*nvec(i:j));
    fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/snps);
  end

  % ensure that freqvec contains frequency of minor allele
  i_major = freqvec > 0.5;
  freqvec(i_major) = 1.0 - freqvec(i_major);

  fname = sprintf('%s/%s_zmat.mat', result_dir, out);
  fprintf('saving %s as -v7.3... ', fname);
  save(fname, '-v7.3', 'zmat_orig', 'zmat_perm', 'beta_orig', 'beta_perm', 'measures', 'nvec', 'zvec_cca', 'freqvec', 'ymat_corr');
  fprintf('OK.\n')
else
  fprintf('loading %s... ', zmat_name);
  load(zmat_name);
  fprintf('OK.\n')
  snps=size(zmat_orig, 1);
  npheno=size(zmat_orig, 2);
end

gwas_time_sec = toc; tic

fprintf('running MOSTest analysis...')
ivec_snp_good = all(isfinite(zmat_orig) & isfinite(zmat_perm), 2);
ivec_snp_good = ivec_snp_good & (freqvec > maf_threshold); % ignore all SNPs with maf < maf_threshold

if use_pheno_corr
  % use correlation structure of the phenotypes
  C0 = ymat_corr;
  C1 = ymat_corr;
else
  % use correlation structure of the z scores, calculated under permutation
  % we don't weight SNPs by LD because the permutation scheme breaks the LD structure
  snps_weight_values = ones(size(zmat_perm, 1), 1);

  % correlation structure of the null z scores
  C0 = weightedcorrs(zmat_perm(ivec_snp_good, :), snps_weight_values(ivec_snp_good));

  % correlation structure of the real z scores
  C1 = weightedcorrs(zmat_orig(ivec_snp_good, :), snps_weight_values(ivec_snp_good)); % & Hvec>0.1 & CRvec>0.95 & max(abs(zmat(:,:,1)),[],1)>abs(norminv(1e-5))),1)');
end

[U S]  = svd(C0); s = diag(S);

%  C0_reg = diag(diag(C0)); % Complete regularization -- results in imperfect gamma fit
%  C0_reg = eye(size(C0)); % Complete regularization -- results in imperfect gamma fit
%  max_lambda = s(min(10, length(s)));
%  max_lambda = min(0.1, s(min(10, length(s)))); % oleksanf: don't regularize unless it's too bad

if (num_eigval_to_keep > 0), max_lambda=s(num_eigval_to_keep); else max_lambda = 0; end;
C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit

%  C0_reg = U*diag(max(s(40),s))*U';
%C0_reg = C0;  % no regularization

mostvecs = NaN(2,snps); minpvecs = NaN(2,snps); maxlogpvecs = NaN(2,snps);
for i  = 1:2
  if i==1, zmat=zmat_orig; else zmat=zmat_perm; end;
  mostvecs(i,:) = dot(inv(C0_reg)*zmat', zmat');
  minpvecs(i,:) = 2*normcdf(-max(abs(zmat), [], 2));
  maxlogpvecs(i, :) = -log10(minpvecs(i, :));
end

[hc_maxlogpvecs hv_maxlogpvecs] = hist(maxlogpvecs(2,ivec_snp_good),1000); chc_maxlogpvecs = cumsum(hc_maxlogpvecs)/sum(hc_maxlogpvecs);
[hc_mostvecs hv_mostvecs] = hist(mostvecs(2,ivec_snp_good),1000); chc_mostvecs = cumsum(hc_mostvecs)/sum(hc_mostvecs);

if use_paretotails
  pd_maxlogpvecs = paretotails(maxlogpvecs(2,ivec_snp_good), 0.0, paretotails_quantile);
  pd_minpvecs_params = upperparams(pd_maxlogpvecs);
  cdf_minpvecs = 1.0 - fixed_paretotails_cdf(pd_maxlogpvecs,hv_maxlogpvecs);
  maxlogpvecs_corr = -log10(fixed_paretotails_cdf(pd_maxlogpvecs, maxlogpvecs));

  pd_mostvecs = paretotails(mostvecs(2,ivec_snp_good),  0.0, paretotails_quantile);
  pd_mostvecs_params = upperparams(pd_mostvecs);
else
  pd_minpvecs = fitdist(colvec(minpvecs(2,ivec_snp_good)),'beta'); % Not a great fit
  pd_minpvecs_params = [pd_minpvecs.a, pd_minpvecs.b];
  cdf_minpvecs=cdf(pd_minpvecs,10.^-hv_maxlogpvecs,'upper');
  maxlogpvecs_corr = -log10(cdf(pd_minpvecs,minpvecs));

  pd_mostvecs = fitdist(colvec(mostvecs(2,ivec_snp_good)),'gamma'); % Seems to work -- beta and wbl  do not
  pd_mostvecs_params = [pd_mostvecs.a, pd_mostvecs.b];
end

if use_paretotails
    cdf_mostvecs = 1.0 - fixed_paretotails_cdf(pd_mostvecs,hv_mostvecs);
    mostvecs_corr = -log10(fixed_paretotails_cdf(pd_mostvecs,mostvecs));
else
    cdf_mostvecs = pd_mostvecs.cdf(hv_mostvecs);
    mostvecs_corr = -log10(cdf(pd_mostvecs,mostvecs,'upper'));
end

fprintf('Done.\n')

fprintf('GWAS yield minP: %d; MOST: %d\n',sum(maxlogpvecs_corr(1,ivec_snp_good)>-log10(5e-8)),sum(mostvecs_corr(1,ivec_snp_good)>-log10(5e-8)));
fprintf('%i\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t\n', npheno, cond(C0), pd_minpvecs_params(1), pd_minpvecs_params(2), pd_mostvecs_params(1), pd_mostvecs_params(2)) 

most_time_sec = toc;

minp_log10pval_orig = maxlogpvecs_corr(1, :);
most_log10pval_orig = mostvecs_corr(1, :);
minp_log10pval_perm = maxlogpvecs_corr(2, :);
most_log10pval_perm = mostvecs_corr(2, :);

% fix potential log10(0) = Inf issues for valid SNPs
ivec_snp_good_flat = reshape(ivec_snp_good.', 1, []);
minp_log10pval_orig(isinf(minp_log10pval_orig) & ivec_snp_good_flat) = -log10(eps(0));
minp_log10pval_perm(isinf(minp_log10pval_perm) & ivec_snp_good_flat) = -log10(eps(0));
most_log10pval_orig(isinf(most_log10pval_orig) & ivec_snp_good_flat) = -log10(eps(0));
most_log10pval_perm(isinf(most_log10pval_perm) & ivec_snp_good_flat) = -log10(eps(0));

fname=sprintf('%s/%s.mat', result_dir, out);
fprintf('saving %s as -v7... ', fname);
save(fname, '-v7', ...
 'most_log10pval_orig', 'minp_log10pval_orig', ...
 'most_log10pval_perm', 'minp_log10pval_perm', ...
 'nvec', 'freqvec', 'ivec_snp_good', ...
 'measures', 'ymat_corr', 'C0', 'C1', ...
 'minpvecs', 'mostvecs', ...
 'hv_maxlogpvecs', 'hc_maxlogpvecs', 'chc_maxlogpvecs', 'cdf_minpvecs', ...
 'hv_mostvecs', 'hc_mostvecs', 'chc_mostvecs', 'cdf_mostvecs', ...
 'pd_minpvecs_params', 'pd_mostvecs_params', 'gwas_time_sec', 'most_time_sec');
fprintf('Done.\n')

fprintf('MOSTest analysis is completed.\n')



if 0
  % QQ plots for minP and MOSTest - permuted vs original distribution of the test statistics (no beta/gamma fitting)
  x=sort(-log10(minpvecs(2, ivec_snp_good')));  y=sort(-log10(minpvecs(1, ivec_snp_good')));  lim=20; % x=permuted; y=original
  figure(1); clf; hold on; plot(x,y, '.'); xlim([0, lim]); ylim([0, lim]); plot([0, lim], [0, lim]); xlabel('minP statistic, permuted'); ylabel('minP statistic, original'); title('minP: QQ plot original vs permuted');
  x=sort(mostvecs(2, ivec_snp_good'));  y=sort(mostvecs(1, ivec_snp_good'));  lim=400;; % x=permuted; y=original
  figure(2); clf; hold on; plot(x,y, '.'); xlim([0, lim]); ylim([0, lim]); plot([0, lim], [0, lim]); xlabel('MOSTest statistic, permuted'); ylabel('MOSTest statistic, original'); title('MOSTest: QQ plot original vs permuted');
end
end % function

