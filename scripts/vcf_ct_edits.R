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
  parser$add_argument('-v', '--vcf', type='character', help='file name of vcf file to be annotated')
  parser$add_argument('-o', '--output', type='character', help='file name of output vcf file')
  args=parser$parse_args()
  
  vcf.filenames = args$vcf
  # only write gz vcf for some reson.....
  output.filename = paste0(args$output,'.gz')
  # vcf.filenames <- '/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/manta_021919/vcf_inv_corrected_dir/C-M916LH-L001-dsomaticSV.vcf'
  # vcf.filenames <- '/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/manta_021919/vcf_inv_corrected_dir/C-001440-L001-dsomaticSV.vcf'
  # output.filename <- '/ifs/work/bergerm1/zhengy1/RET_all/Code/tmp.vcf.gz'
  vcf.file <- read.vcfR(vcf.filenames,verbose = F)
  vcf.data <- data.table(vcf.file@fix)
  if(nrow(vcf.data) > 0){
    # mapping of read directions to alt annotation
    CT.vector <- structure(c('5to3','3to3','3to5','5to5'),names = c('t[p[','t]p]',']p]t','[p[t'))
    vcf.data[,c('EventType','BND_CT') := list(gsub('Manta','',str_extract(ID,'Manta...')),
                                              ifelse(grepl('*.\\[.*.\\[',ALT),'t[p[',
                                                     ifelse(grepl('*.\\].*.\\]',ALT),'t]p]',
                                                            ifelse(grepl('\\].*.\\].*',ALT),']p]t',
                                                                   ifelse(grepl('\\[.*.\\[.*',ALT),'[p[t',NA)))))] %>% rowwise %>%
      # delection insertion duplication always 5 to 3
      mutate(CT = ifelse(EventType %in% c('DEL','INS','DUP'),'5to3',
                         # inversion either 3to3 or 5to5 (use INV3/5 tag in INFO)
                         ifelse(EventType == 'INV',ifelse(grepl('INV3',INFO),'3to3','5to5'),
                                # mapping of translocation event type
                                ifelse(EventType == 'BND',CT.vector[[BND_CT]],'NA'))),
             # extract mate chr and end pos
             mate.ID = gsub('MATEID=|;','',str_extract(INFO,'MATEID=Manta...:[0-9:]+;'))) %>%
      merge(vcf.data[,.(CHROM,POS,ID)],by.x = 'mate.ID',by.y = 'ID',all.x = T,suffixes = c('','.mate'))  %>%
      mutate(END = ifelse(is.na(POS.mate),'',paste0('END=',POS.mate)),CHR2 = paste0('CHR2=',ifelse(is.na(CHROM.mate),CHROM,CHROM.mate))) %>%
      mutate(INFO = paste0(INFO,';CT=',CT,';',CHR2,';',END)) %>% data.table() -> vcf.data
    # rows with GLxxxxx as chromosome
    row.to.del <- c(grep('GL',vcf.data$CHROM),grep('GL',vcf.data$CHROM.mate))
    if(length(row.to.del) > 0){
      vcf.file@fix <- as.matrix(vcf.data[-row.to.del,])
      vcf.file@gt <- as.matrix(data.table(vcf.file@gt)[-row.to.del,])
    }else{
      vcf.file@fix <- as.matrix(vcf.data)
    }
    vcf.data <- vcf.data %>% select(-c(EventType,CT,BND_CT,CHROM.mate,mate.ID,POS.mate,END,CHR2)) %>% data.table()
    
    if(any(grepl('CT=NA',vcf.data$INFO))){
      warning(paste0('there are unknown connection type in this vcf file -- ',vcf.filenames))
    }
  }
  write.vcf(x = vcf.file,file = output.filename)
  # process and discard gz file
  system(paste0('zcat ',output.filename,' > ',gsub('.gz$','',output.filename)))
  unlink(output.filename)
}
