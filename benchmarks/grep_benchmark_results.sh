experiment_dir=$1
HEAD="algorithm,graph,timeout,seed,k,epsilon,num_threads,imbalance,totalPartitionTime,objective,km1,cut,failed"

mkdir experimental_results
for result_folder in $experiment_dir/*_results;
do
  result_file="${result_folder/_results/.csv}"
  header_file="${result_folder/_results/.header.csv}"
  rm -f $result_file
  if [[ -f "$header_file" ]]; then
    cat "$header_file" >> $result_file
  else
    echo "$HEAD" >> $result_file
  fi
  for instance_result_file in $result_folder/*.results;
  do
    tail -1 $instance_result_file >> $result_file
  done
  mv $result_file experimental_results/
done

tar -cvf experimental_results.tar experimental_results/*
rm -rf experimental_results