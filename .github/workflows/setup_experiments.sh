declare -a blocks=("2" "8")

for seed in $(seq 3)
do
  for instance in benchmarks/*.hgr;
  do
    instance_name=`echo $instance | sed 's!.*/!!'`
    for k in "${blocks[@]}"
    do
      echo "\$1 -h $instance -k $k -e 0.03 -o km1 -m direct -p \$2 --seed=$seed --sp-process=false --csv=true --verbose=false --enable-progress-bar=false --show-detailed-timings=true --s-num-threads=2 >> \$3/$instance_name.$k.2.$seed.results" >> benchmark.sh
    done
  done
done
chmod +x benchmark.sh
