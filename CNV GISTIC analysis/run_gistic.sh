#!/bin/sh
## run example GISTIC analysis

## output directory
echo --- creating output directory ---
basedir=`pwd`/result
mkdir -p $basedir 

echo --- running GISTIC ---
## input file definitions
segfile=`pwd`/CNV_155_gistic_res_change.txt
#markersfile=`pwd`/examplefiles/markersfile.txt
refgenefile=`pwd`/hg19.mat
#alf=`pwd`/examplefiles/arraylistfile.txt
#cnvfile=`pwd`/examplefiles/cnvfile.txt
## call script that sets MCR environment and calls GISTIC executable 
`pwd`/bin/GISTIC/gistic2 -b $basedir -seg $segfile -refgene $refgenefile -ta 0.585 -td -1 -broad 1 -brlen 0.5 -maxseg 10000 -conf 0.90 -genegistic 1 -armpeel 1 -smallmem 1 -savegene 1 
