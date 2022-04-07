#!/bin/bash

module load gcc/6.2.0
module load stacks/2.3e

#gstacks
gstacks -I $src/data/ -M $src/STACKS/popmap -O $src/STACKS

#populations
populations -P $src/STACKS -M $src/STACKS/pop_with_order.txt -O $src/STACKS/results -t 8 \
--min-mac 1 -R 0.2 --write-single-snp \
-r 0.9 --min-maf 0.05 --fstats  --vcf
#--genepop