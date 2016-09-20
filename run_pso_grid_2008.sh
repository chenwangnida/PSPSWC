 #!/bin/sh

need sgegrid

NUM_RUNS=50

for i in {1..1}; do
  qsub -t 1-$NUM_RUNS:1 graph_pso.sh ~/workspace/PSPSWC/WSC08TestSet0${i} WSC08SWSC${i};
done