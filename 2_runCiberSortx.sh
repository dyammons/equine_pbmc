#!/usr/bin/bash

#Run with no batch correction
singularity exec \
-B /pl/active/dow_lab/dylan/eq_long_scRNA/pbmc/output/cibersort_data:/src/data \
-B /pl/active/dow_lab/dylan/eq_long_scRNA/pbmc/output/cibersort_output/noBatch:/src/outdir \
/pl/active/dow_lab/dylan/software/cibersort/fractions_latest.sif \
/src/CIBERSORTxFractions \
--username dyammons@colostate.edu \
--token b451a793ee4a3357367ea2fb3d551fcc \
--single_cell TRUE \
--refsample scrna_celltype_cibersort.txt \
--mixture mixture.txt \
--fraction 0 \
--replicates 30 \
--remake TRUE

#Run with bulk mode batch correction
singularity exec \
-B /pl/active/dow_lab/dylan/eq_long_scRNA/pbmc/output/cibersort_data:/src/data \
-B /pl/active/dow_lab/dylan/eq_long_scRNA/pbmc/output/cibersort_output/bBatch:/src/outdir \
/pl/active/dow_lab/dylan/software/cibersort/fractions_latest.sif \
/src/CIBERSORTxFractions \
--username dyammons@colostate.edu \
--token b451a793ee4a3357367ea2fb3d551fcc \
--single_cell TRUE \
--refsample scrna_celltype_cibersort.txt \
--mixture mixture.txt \
--fraction 0 \
--replicates 30 \
--remake TRUE \
--rmbatchBmode TRUE

#Run with single cell batch correction, the reccomened correction method
singularity exec \
-B /pl/active/dow_lab/dylan/eq_long_scRNA/pbmc/output/cibersort_data:/src/data \
-B /pl/active/dow_lab/dylan/eq_long_scRNA/pbmc/output/cibersort_output/sBatch:/src/outdir \
/pl/active/dow_lab/dylan/software/cibersort/fractions_latest.sif \
/src/CIBERSORTxFractions \
--username dyammons@colostate.edu \
--token b451a793ee4a3357367ea2fb3d551fcc \
--single_cell TRUE \
--refsample scrna_celltype_cibersort.txt \
--mixture mixture.txt \
--fraction 0 \
--replicates 30 \
--remake TRUE \
--rmbatchSmode TRUE
