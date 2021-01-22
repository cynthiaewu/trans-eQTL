#!/bin/bash

OUTDIR=/storage/mgymrek/cpma/power/
SSIZE=500
TARGETS="5 15 20 30 40 60 80 100 150 200 250 300 350 400 450 500 700 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000"
BETAS="001 002 003 004 005 01 02 03 04 05 1"
RESDIR=/gymreklab-tscc/cynthiawu/Test_nullsnps/simulate_eqtls_only/FastMultivariate/Single_eqtl

for targets in $TARGETS
do
    for beta in $BETAS
    do
        # CPMA
        echo $targets $beta
#        cat ${RESDIR}/SampleSize${SSIZE}/SingleParameter/numTarget_${targets}/Beta_${beta}/Simulation_*/CPMAx/gene-snp-eqtl_cpmax_pvalues_1.0 | \
#            grep SNP0 | cut -f 4 > ${OUTDIR}/CPMA_${SSIZE}_${targets}_${beta}.pvals
        
        # Mixture model
        cat ${RESDIR}/SampleSize${SSIZE}/SingleParameter/numTarget_${targets}/Beta_${beta}/Simulation_*/mixtureModel/gene-snp-eqtl_mixturepvalue | \
            grep SNP0 | cut -f 2 > ${OUTDIR}/Mixture_${SSIZE}_${targets}_${beta}.pvals
        
        # Pairwise
#        for f in $(ls ${RESDIR}/SampleSize${SSIZE}/SingleParameter/numTarget_${targets}/Beta_${beta}/Simulation_*/CPMA/gene-snp-eqtl)
#	do 
#	    cat $f | grep -v gene | datamash min 5
#	done > ${OUTDIR}/Pairwise_${SSIZE}_${targets}_${beta}.pvals

    done
done
