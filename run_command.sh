#!/bin/bash

./nextflow run main.nf \
-c resources.config \
--obj /home/lgdqamorim/scratch/projeto_nextflow/data/vehicle_final_server.rds \
--column clone_annotation \
--species human \
--n_cells 3000 \
--n_cores 32 \
--target /home/lgdqamorim/scratch/projeto_nextflow/alvos.txt \
--network hdwgcna \
--outdir results \
-profile singularity -resume
