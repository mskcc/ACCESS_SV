process.support <- function(evidence.type,evidence.data){
  # evidence.type <- vcf.genotype$FORMAT[1]
  # evidence.data <- vcf.genotype$value[1]
  paste0(paste0(unlist(strsplit(evidence.type,':')[[1]]),'=',unlist(strsplit(evidence.data,':')[[1]])),collapse = '--')
}



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
  parser$add_argument('-t', '--tab', type='character', help='file name of tab file to be annotated')
  parser$add_argument('-v', '--vcf', type='character', help='file name of vcf file to be annotated')
  parser$add_argument('-o', '--output', type='character', help='file name of output file')
  args=parser$parse_args()
  
  vcf.filenames = args$vcf
  tab.filenames = args$tab
  output.filename = args$output
  prefix = args$prefix
  # vcf.filenames <- '/ifs/work/bergerm1/RNAseq/Analysis/manta/run_030619/vcf_dir/Sample_P-0013196-T04-RNAS_IGO_05500_FO_1rnaSV.vcf'
  vcf.filenames <- '/ifs/work/bergerm1/zhengy1/RET_all/Code/test/SV_Test/somaticSV_inv_corrected_edited.vcf'
  tab.filenames = '/ifs/work/bergerm1/zhengy1/RET_all/Code/test/SV_Test/C-02WK6K-L001-d_Annotated.txt'
  output.filename <- '/ifs/work/bergerm1/zhengy1/RET_all/Code/test/SV_Test/C-02WK6K-L001-d_Annotated_Evidence.txt'
  tab.data <- fread(tab.filenames) %>% mutate(chr1 = as.character(chr1),pos1 = as.numeric(pos1),chr2 = as.character(chr2),pos2 = as.numeric(pos2)) %>% data.table()
  vcf.file <- read.vcfR(vcf.filenames,verbose = F)
  vcf.data <- data.table(vcf.file@fix)
  gene.list <- c('ALK','BRAF','ERBB2','EGFR','FGFR1','FGFR2','FGFR3','KIT','MET','NTRK1','NTRK2','NTRK3','PDGFRA','RET','ROS1')
  if(nrow(vcf.data)== 0){
    output.data <- data.frame(matrix(nrow = 0,ncol = 37))
    colnames(output.data) = c('chr1','pos1','chr2','pos2','Connection_Type','SV_LENGTH','SV_Type','BrkptType','PairedEndSupport.nomal',
                              'PairedEndRef.nomal','SplitReadSupport.nomal','SplitReadRef.nomal','PairedEndSupport.tumor',
                              'PairedEndRef.tumor','SplitReadSupport.tumor','SplitReadRef.tumor','str1','str2','gene1',
                              'transcript1','site1','gene2','transcript2','site2','fusion','Cosmic_Fusion_Counts',
                              'repName-repClass-repFamily:-site1','repName-repClass-repFamily:-site2','CC_Chr_Band',
                              'CC_Tumour_Types(Somatic)','CC_Cancer_Syndrome','CC_Mutation_Type','CC_Translocation_Partner',
                              'DGv_Name-DGv_VarType-site1','DGv_Name-DGv_VarType-site2','TumorId','NormalId')
  }else{
    # rows with GLxxxxx as chromosome
    row.to.del <- grep('GL',vcf.data$CHROM)
    if(length(row.to.del) > 0){
      vcf.file@fix <- as.matrix(vcf.data[-row.to.del,])
      vcf.file@gt <- as.matrix(data.table(vcf.file@gt)[-row.to.del,])
    }
    vcf.genotype <- data.table(vcf.file@gt) %>% cbind(vcf.data[,.(chr1 = CHROM,pos1 = as.numeric(POS),
                                                                  chr2 = gsub('CHR2=|;','',str_extract(INFO,'CHR2=[0-9XYM]+;')),
                                                                  pos2 = as.numeric(gsub('END=|;','',str_extract(INFO,'END=[0-9]+;|END=[0-9]+$'))),
                                                                  Connection_Type = gsub('CT=|;','',str_extract(INFO,'CT=[0-9]to[0-9];')),
                                                                  SV_LENGTH = ifelse(grepl('SVLEN',INFO),
                                                                                     as.numeric(gsub('SVLEN=|;','',str_extract(INFO,'SVLEN=[0-9]+;|SVLEN=-[0-9]+;'))),0),
                                                                  SV_Type = gsub('SVTYPE=|;','',str_extract(INFO,'SVTYPE=[A-Z]+;')),INFO
                                                                  )]) %>% 
      mutate(BrkptType = ifelse(grepl('IMPRECISE',INFO),'IMPRECISE',ifelse(grepl('SR',FORMAT),'PRECISE','IMPRECISE'))) %>% select(-INFO) %>% data.table()
    normal.ID <- colnames(vcf.genotype)[2]
    tumor.ID <- colnames(vcf.genotype)[3]
    vcf.genotype %>% melt(id.vars = c('chr1','pos1','chr2','pos2','FORMAT','Connection_Type','SV_LENGTH','SV_Type','BrkptType')) %>% rowwise() %>% 
      mutate(Evidence = process.support(FORMAT,value)) %>%
      mutate(sample.type = ifelse(variable == normal.ID,'Normal','Tumor'),
             PairedEndSupport = ifelse(grepl('PR',Evidence),strsplit(gsub('PR=','',str_extract(Evidence,'PR=[0-9]+,[0-9]+')),',')[[1]][2]),
             PairedEndRef = ifelse(grepl('PR',Evidence),strsplit(gsub('PR=','',str_extract(Evidence,'PR=[0-9]+,[0-9]+')),',')[[1]][1]),
             SplitReadSupport = ifelse(grepl('PR',Evidence),strsplit(gsub('SR=','',str_extract(Evidence,'SR=[0-9]+,[0-9]+')),',')[[1]][2]),
             SplitReadRef = ifelse(grepl('PR',Evidence),strsplit(gsub('SR=','',str_extract(Evidence,'SR=[0-9]+,[0-9]+')),',')[[1]][1])) %>%
      select(-c(FORMAT,variable,value,Evidence)) %>% data.table() -> vcf.genotype
    output.data <- merge(vcf.genotype[sample.type == 'Normal',!'sample.type',with = F],vcf.genotype[sample.type == 'Tumor',!'sample.type',with = F],
                         by = c('chr1','pos1','chr2','pos2','Connection_Type','SV_LENGTH','SV_Type','BrkptType'),all = T,suffixes = c('.nomal','.tumor')) %>%
      merge(tab.data,by = c('chr1','pos1','chr2','pos2'),all = T) %>% rowwise %>% 
      mutate(TumorId = tumor.ID,NormalId = normal.ID,`Sign out` = ifelse(gene1 %in% gene.list | gene2 %in% gene.list,'Y','N')) %>%
      select(TumorId,NormalId,`Sign out`,everything()) %>%
      data.table()
  }
  write.table(output.data,output.filename,sep = '\t',quote = F,row.names = F)
  print(paste0('Writing to -- ',output.filename))
}
