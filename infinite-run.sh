#!/bin/bash

while true; do 
  TIMESTAMP=$(date +"%T.%N")
  OUTPUT="$(./build/mt-kahypar/application/MtKaHyPar.exe -h ./tests/instances/ibm01.hgr --preset-type=default -k 8 -e 0.03 -t 10 -m direct -okm1 --seed=1 | grep -E '^.*km1\s+= [0-9]+ \(primary objective function\)' | awk '{print $3;}')"
  echo "$TIMESTAMP, $OUTPUT" >> output.csv
done
