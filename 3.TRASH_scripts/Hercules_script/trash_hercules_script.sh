#!/bin/bash
#$ -cwd
#$ -V
#$ -N 'trash'
#$ -q h0809.q


source ~/.bash_profile
conda activate TRASH_env


export PATH=/users-d3/vadim.pisarenco/Programs_V/TRASH_gitcloned/:$PATH
which TRASH_run.sh

date
GENOME=/users-d3/vadim.pisarenco/Programs_V/TRASH_gitcloned/example_run/CP068268_39050443_39150442.fa #ABSOLUTE PATH !!!!!!!!!
TRASH_run.sh $GENOME --def
date

