#!/bin/sh
#
# Force Bourne Shell if not Sun Grid Engine default shell (you never know!)
#
#$ -S /bin/sh
#
# I know I have a directory here so I'll use it as my initial working directory
#
#$ -wd /vol/grid-solar/sgeusers/wangchen
#
# End of the setup directives
#
# Stdout from programs and shell echos will go into the file
#    scriptname.o$JOB_ID
#  so we'll put a few things in there to help us see what went on
#

# For testing locally
#JOB_ID=3
#SGE_TASK_ID=1

DIR_TMP="/local/tmp/wangchen/$JOB_ID/"
#DIR_HOME="/am/courtenay/home1/wangchen/"
DIR_HOME="/home/wangchen/"
#DIR_GRID=$DIR_HOME"grid/"
DIR_GRID="/vol/grid-solar/sgeusers/wangchen/"
DIR_WORKSPACE="workspace/"
DIR_PROGRAM=$DIR_HOME$DIR_WORKSPACE"PSPSWC/"
DIR_OUTPUT=$DIR_GRID$2 # Match this argument with dataset name

FILE_JOB_LIST="CURRENT_JOBS.txt"
FILE_RESULT_PREFIX="out"

mkdir -p $DIR_TMP

# copy the local lib to grid
#mkdir -p lib
#cp ~/lib/*  .

# Preliminary test to ensure that the directory has been created successfully.
if [ ! -d $DIR_TMP ]; then
  echo "Could not create the temporary directory for processing the job. "
  echo "/local/tmp/ directory: "
  ls -la /local/tmp
  echo "/local/tmp/wangchen directory: "
  ls -la /local/tmp/wangchen
  echo "Exiting"
  exit 1
fi

# Add the job ID to the list of jobs currently being processed.
echo $JOB_ID >> $DIR_GRID$FILE_JOB_LIST

# Copy the files required for processing into the temporary directory.
cp -r $DIR_PROGRAM"bin" $DIR_TMP
cp $1/* $DIR_TMP # Copy datasets
cp ~/lib/* $DIR_TMP # Copy jars

mkdir -p $DIR_TMP"results"

cd $DIR_TMP

# Give me permission to read, write, and execute everything in my temp directory
chmod +700 *

echo "Running: "

seed=$SGE_TASK_ID
result=$FILE_RESULT_PREFIX$seed.stat

java -classpath ./bin:.:jgraph-5.13.0.0.jar:jgrapht-core-0.9.2.jar:jgrapht-ext-0.9.2-uber.jar:jgrapht-ext-0.9.2.jar:jgraphx-2.0.0.1.jar ec.pso.GraphPSO $result problem.xml services-output.xml taxonomy.owl $seed

cp $result ./results

# Now we move the output to a place to pick it up from later and clean up
cd results
if [ ! -d $DIR_OUTPUT ]; then
  mkdir $DIR_OUTPUT
fi
cp $DIR_TMP"results"/*.stat $DIR_OUTPUT

# Do the cleaning up from our starting directory
rm -rf /local/tmp/wangchen/$JOB_ID

# Remove the job ID from the list of jobs currently being processed.
FILE_TMP=`date +%N`
grep -v "$JOB_ID" $DIR_GRID$FILE_JOB_LIST > $DIR_GRID$FILE_TMP
mv $DIR_GRID$FILE_TMP $DIR_GRID$FILE_JOB_LIST

# Finish the job.
echo "Ran through OK"


