experiment_dir=$1
HEAD="algorithm,threads,graph,k,seed,epsilon,imbalance,objective,km1,cut,partitionTime,fmTime,lpTime,coarseningTime,ipTime,preprocessingTime"

for result_folder in $experiment_dir/*/;
do
  folder_name=`basename "$result_folder"`
  echo "$HEAD" >> $experiment_dir/$folder_name.csv
  for instance_result_file in $result_folder/*.results;
  do
    tail -1 $instance_result_file >> $experiment_dir/$folder_name.csv
  done
done