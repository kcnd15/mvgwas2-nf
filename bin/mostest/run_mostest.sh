# run MOSTest-script with MATLAB
# kc
# arguments: pheno_file, bfile_prefix, out_prefix, data_dir, result_dir
# ./run_mostest.sh pheno.txt chr21 mostest_results ../../data/mostest ../../result

echo "run_mostest.sh command line arguments:"
echo "try, run('run_mostest('$1', '$2', '$3', '$4', '$5')'; end; quit"
pwd
# cd /home/kcan/UOC/tfm/mvgwas2-nf/bin/mostest

export MATLABPATH=/home/kcan/UOC/tfm/mvgwas2-nf/bin/mostest

# execute matlab script run_mostest.m
matlab -nosplash -nodesktop -r "try, run_mostest('$1', '$2', '$3', '$4', '$5'), catch me, fprintf('run_mostest.sh: %s: %s\n', me.identifier, me.message), exit(1);, end, exit(0);"
