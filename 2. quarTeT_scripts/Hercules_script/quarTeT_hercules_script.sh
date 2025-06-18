#!/bin/bash
#$ -cwd
#$ -V
#$ -N quartet_script
#$ -q h0107.q
#$ -pe ompi8 1

source ~/.bashrc

conda activate quartet

genome=/users-d2/guillem.garcia/genomes/rhynchospora_pubera/ncbi_dataset/data/GCA_028095005.1/GCA_028095005.1_MPIPZ_rhyPub2m_2.0_genomic.fna

q_path=/users-d2/guillem.garcia/quarTeT/quartet.py

python3 $q_path CentroMiner -i $genome
