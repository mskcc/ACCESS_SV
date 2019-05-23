#!/bin/bash

bam_file=$1
#bam_file=/ifs/work/bergerm1/brannona/ACCESS_M1.8/ACCESSv1-VAL-20180034/standard/F19_cl_aln_srt_MD_IR_FX_BR.bam
output_file=$2

samtools view -h $bam_file \
	| awk '$6!~/^.*.[0-9]+I[0-9]+D$|^.*.[0-9]+D[0-9]+I$|^.*.[0-9]+D[0-9]+I[0-9]+S$|^.*.[0-9]+I[0-9]+D[0-9]+S$/' \
	| samtools view -Sb - > $output_file
