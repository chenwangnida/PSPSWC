 #!/bin/sh

need sgegrid

NUM_RUNS=40

for i in {1..5}; do
  qsub -t 1-$NUM_RUNS:1 graph_pso.sh ~/workspace/swsc2009/Set0${i}MetaData 2009-semantic-QoS-warepso${i};
done