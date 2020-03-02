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
  # this library writes objects of class vcfR to *.vcf.gz format only.
  output.filename = args$output
  output.filename.gzipped = paste0(output.filename, '.gz')
  vcf.file <- read.vcfR(vcf.filenames,verbose = F)
  vcf.data <- data.table(vcf.file@fix)
  # if the input vcf has one or more variants
  if(nrow(vcf.data) > 0){
    # mapping of read directions to alt annotation
    CT.vector <- structure(c('5to3','3to3','3to5','5to5'),names = c('t[p[','t]p]',']p]t','[p[t'))
    vcf.data[,c('EventType','BND_CT') := list(gsub('Manta','',str_extract(ID,'Manta...')),
                                              ifelse(grepl('*.\\[.*.\\[',ALT),'t[p[',
                                                     ifelse(grepl('*.\\].*.\\]',ALT),'t]p]',
                                                            ifelse(grepl('\\].*.\\].*',ALT),']p]t',
                                                                   ifelse(grepl('\\[.*.\\[.*',ALT),'[p[t',NA)))))] %>% rowwise %>%
      # delection insertion duplication always 5 to 3
      mutate(CT = ifelse(EventType %in% c('INS','DUP'),'5to3',
                         ifelse(EventType == 'DEL','3to5',
                                # inversion either 3to3 or 5to5 (use INV3/5 tag in INFO)
                                ifelse(EventType == 'INV',ifelse(grepl('INV3',INFO),'3to3','5to5'),
                                       # mapping of translocation event type
                                       ifelse(EventType == 'BND',CT.vector[[BND_CT]],'NA')))
                         ),
             # extract mate chr and end pos
             mate.ID = gsub('MATEID=|;','',str_extract(INFO,'MATEID=Manta...:[0-9:]+;'))) %>%
      merge(vcf.data[,.(CHROM,POS,ID)],by.x = 'mate.ID',by.y = 'ID',all.x = T,suffixes = c('','.mate'))  %>%
      mutate(END = ifelse(is.na(POS.mate),'',paste0('END=',POS.mate)),CHR2 = paste0('CHR2=',ifelse(is.na(CHROM.mate),CHROM,CHROM.mate))) %>%
      mutate(INFO = paste0(INFO,';CT=',CT,';',CHR2,';',END)) %>% data.table() -> vcf.data
    # rows with GLxxxxx as chromosome
    row.to.del <- c(grep('GL',vcf.data$CHROM),grep('GL',vcf.data$CHROM.mate))
    vcf.data <- vcf.data %>% select(-c(EventType,CT,BND_CT,CHROM.mate,mate.ID,POS.mate,END,CHR2)) %>% data.table()
    
    # if none of the variants belong to autosomes or sex chromosome,
    #  remove the variants and manually write the vcf header to output vcf file.
    #  This is because an empty matrix cannot be assigned to vcf.file@gt at line
    #  `vcf.file@fix <- as.matrix(vcf.data[-row.to.del,])`
    if(length(row.to.del) == nrow(vcf.data)) {
	    temp_vcf <- readLines(vcf.filenames)
	    nth.row <- length(vcf.file@meta)+1
	    fileHandle <- file(output.filename)
	    writeLines(temp_vcf[1:nth.row], fileHandle)
	    close(fileHandle)
	    # Quit proccess after manually writing the output vcf file
	    quit(status=0, save="no")
    }
    # if there are variants left after removing non-autosomes and non-sex_chromosome
    #  variants, proceed with using vcfR methods for writing output
    else {
	if(length(row.to.del) > 0){
      	    vcf.file@fix <- as.matrix(vcf.data[-row.to.del,])
	    vcf.file@gt <- as.matrix(data.table(vcf.file@gt)[-row.to.del,])
    	}else{
            vcf.file@fix <- as.matrix(vcf.data)
   	}

    	if(any(grepl('CT=NA',vcf.data$INFO))){
      		warning(paste0('there are unknown connection type in this vcf file -- ',vcf.filenames))
    	}
    }
  }

  # Finally write the output vcf and unzip it
  write.vcf(x = vcf.file, file = output.filename.gzipped)
  # process and discard gz file
  system(paste0('zcat ', output.filename.gzipped, ' > ', output.filename))
  unlink(output.filename.gzipped)
}
