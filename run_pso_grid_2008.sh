 #!/bin/sh

need sgegrid

NUM_RUNS=50

for i in {1..1}; do
  qsub -t 1-$NUM_RUNS:1 graph_pso.sh ~/workspace/swsc2008/Set0${i}MetaData 2008-semantic-pso${i};
done
