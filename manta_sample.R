# Executable -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
  library(stringr)
  library(vcfR)
  library(dplyr)
})

if (!interactive()) {
  
  parser=ArgumentParser()
  parser$add_argument('-t', '--tumor', type='character', help='tumor bam for SV calling')
  parser$add_argument('-n', '--normal', type='character', help='normal bam for SV calling')
  parser$add_argument('-o', '--output', type='character', help='output directory for manta')
  parser$add_argument('-f', '--fasta', type='character', help='fasta file for alignment')
  args=parser$parse_args()
  
  tumor_path = args$tumor
  normal_path = args$normal
  manta.dir = args$output
  fasta = args$fasta
  
  system(paste0(
    '/opt/common/CentOS_6-dev/manta/1.5.0/bin/configManta.py --normalBam ',
    normal_path ,' --tumorBam ',tumor_path,
    ' --runDir ',manta.dir,' --exome ',
    '--referenceFasta ',fasta
  ))
  system(paste0(
    'python ',manta.dir,'/runWorkflow.py -m local -j 8'
  ))
}
