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
  parser$add_argument('-t', '--tumor', type='character', help='Tumor bam for SV calling')
  parser$add_argument('-n', '--normal', type='character', help='Normal bam for SV calling')
  parser$add_argument('-o', '--output', type='character', help='Output directory for manta')
  parser$add_argument('-f', '--fasta', type='character', help='Fasta file for alignment')
  parser$add_argument('-m', '--manta', type='character', help='Directory of manta installation')
  args=parser$parse_args()
  
  tumor_path = args$tumor
  normal_path = args$normal
  out.dir = args$output
  fasta = args$fasta
  manta.dir = args$manta
  
  system(paste0(
    manta.dir,'/bin/configManta.py --normalBam ',
    normal_path ,' --tumorBam ',tumor_path,
    ' --runDir ',out.dir,' --exome ',
    '--referenceFasta ',fasta
  ))
  system(paste0(
    'python ',out.dir,'/runWorkflow.py -m local -j 8'
  ))
}
