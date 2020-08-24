#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
})

if (!interactive()) {
  parser=ArgumentParser()
  parser$add_argument('-v', '--vcf', type='character', help='file name of vcf file to be annotated')
  parser$add_argument('-o', '--output', type='character', help='file name of output vcf file')
  
  args=parser$parse_args()
  vcf = data.table::fread(args$vcf, colClasses = "character", 
                          sep="\t", header = F)
  vcf_meta = vcf[grepl("^#", vcf$V1)]$V1
  vcf_main = vcf[!grepl("^#", vcf$V1)]$V1
  
  # write vcf metadata
  write.table(vcf_meta, args$output, sep = "\t", col.names = F,
              row.names = F, quote = F, append = F)
  
  if (length(vcf_main) > 0) {
    vcf_data = as.data.table(tstrsplit(vcf_main, "\t", fixed=T))
    names(vcf_data) = strsplit(last(vcf_meta), "\t", fixed = T)[[1]]
    
    # vector mapping 5'/3' connection type to mantas output
    CT.vector = structure(c('5to3','3to3','3to5','5to5'),
                          names = c('t[p[','t]p]',']p]t','[p[t'))
    
    # Generate eventype and bnd_ct columns
    vcf_data[, c('EventType','BND_CT') := list(
      #gsub('Manta','',str_extract(ID,'Manta...')),
      gsub('Manta', '', sub("(Manta...).*", "\\1", ID)),
      ifelse(grepl('*.\\[.*.\\[', ALT),'t[p[',
             ifelse(grepl('*.\\].*.\\]', ALT),'t]p]',
                    ifelse(grepl('\\].*.\\].*', ALT),']p]t',
                           ifelse(grepl('\\[.*.\\[.*', ALT),'[p[t', 'NA'))))),
      by = seq_len(NROW(vcf_data))]
    
    # infer 5'/3' connection type and mate.ID
    vcf_data[, c('CT', 'mate.ID') := list(
      ifelse(EventType %in% c('INS','DUP'), '5to3', 
             ifelse(EventType == 'DEL', '3to5', 
                    # inversion either 3to3 or 5to5 (use INV3/5 tag in INFO)
                    ifelse(EventType == 'INV', 
                           ifelse(grepl('INV3',INFO), '3to3', '5to5'),
                           # mapping of translocation event type
                           ifelse(EventType == 'BND', CT.vector[[BND_CT]], 'NA')))),
      gsub('MATEID=|;', '', sub(".*;(MATEID=Manta...:[0-9:]+);.*", "\\1", INFO))), 
      by = seq_len(NROW(vcf_data))]
    
    # merge data.table with itself on ID and mate.ID to match
    #  translocation partners
    vcf_data = merge(vcf_data, vcf_data, 
                     by.x="ID", by.y = "mate.ID", all.x=T, all.y=F, 
                     suffixes = c('', '.mate'))
    
    vcf_data[, c("END", "CHR2") := list(
      ifelse(is.na(POS.mate), '', paste0('END=',POS.mate)),
      paste0('CHR2=', ifelse(is.na(`#CHROM.mate`), `#CHROM`, `#CHROM.mate`))),
      by = seq_len(NROW(vcf_data))]
    
    vcf_data[, INFO := paste0(INFO,';CT=',CT,';',CHR2,';',END), by = seq_len(NROW(vcf_data))]
    
    # Select for valid Chromosomes
    Chromosomes = c(as.character(1:22), "X", "Y", "MT")
    vcf_data = vcf_data[which(vcf_data$`#CHROM` %in% Chromosomes),]
    
    # Select vcf columns
    genotype_columns_to_keep = grep(
      paste(c(".mate$", "mate.ID", "EventType", 
              "BND_CT", "CT", "END", "CHR2", "ID", 
              "#CHROM", "POS", "REF", "ALT", "QUAL", 
              "FILTER", "INFO", "FORMAT"), collapse = "|"), 
      names(vcf_data), value = T, invert = T)
    columns_to_keep = c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
      "FILTER", "INFO", "FORMAT", genotype_columns_to_keep)
    vcf_data = vcf_data[, ..columns_to_keep]
    
    # Reorder rows
    vcf_data$`#CHROM` = factor(vcf_data$`#CHROM`, levels=Chromosomes)
    vcf_data = vcf_data[order(`#CHROM`, as.numeric(POS))]
    
    # Warn about unknown connection types
    if (any(grepl('CT=NA',vcf_data$INFO))) {
      warning(paste0(
        'There are unknown connection type in this vcf file -- ', args.vcf))
    }
    
    # Write vcf data
    write.table(vcf_data, args$output, sep = "\t", row.names = F,
                col.names = F, quote = F, append = T)
  }
}





