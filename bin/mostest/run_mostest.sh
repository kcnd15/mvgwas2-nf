# run MOSTest-script with MATLAB
# kc

echo "run_mostest.sh command line arguments: $1"
echo "try, run('run_mostest('$1', '$2', '$3')'; end; quit"

# matlab -nosplash -nodesktop -r "try, run('run_mostest('$1', '$2', '$3')'; end; quit"
matlab -nosplash -nodesktop -r "try, run_mostest('$1', '$2', '$3'), catch me, fprintf('%s: %s\n', me.identifier, me.message), exit(1);, end, exit(0);"
