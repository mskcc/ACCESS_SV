#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
})


# Constants
sig_gene_list = c(
  'ALK','BRAF','ERBB2','EGFR','FGFR1', 'FGFR2','FGFR3',
  'KIT','MET','NTRK1', 'NTRK2','NTRK3','PDGFRA','RET','ROS1')

access_gene_list = c(
  'AKT1','ALK','APC','AR','ARAF', 'ARID1A','ARID2','ATM','B2M',
  'BCL2','BCOR','BRAF','BRCA1','BRCA2','CARD11','CBFB', 'CCND1',
  'CDH1','CDK4','CDKN2A','CIC','CREBBP','CTCF','CTNNB1','DICER1',
  'DIS3','DNMT3A','EGFR','EIF1AX','EP300','ERBB2', 'ERBB3',
  'ERCC2','ESR1','ETV6','EZH2','FBXW7','FGFR1','FGFR2','FGFR3',
  'FGFR4','FLT3','FOXA1','FOXL2','FOXO1','FUBP1','GATA3', 'GNA11',
  'GNAQ','GNAS','H3F3A','HIST1H3B','HRAS','IDH1','IDH2','IKZF1',
  'INPPL1','JAK1','KDM6A','KEAP1','KIT','KNSTRN','KRAS', 'MAP2K1',
  'MAPK1','MAX','MDM2','MED12','MET','MLH1','MSH2','MSH3','MSH6',
  'MTOR','MYC','MYCN','MYD88','MYOD1','NF1','NFE2L2', 'NOTCH1',
  'NRAS','NTRK1','NTRK2','NTRK3','NUP93','PAK7','PDGFRA','PIK3CA',
  'PIK3CB','PIK3R1','PIK3R2','PMS2','POLE','PPP2R1A', 'PPP6C',
  'PRKCI','PTCH1','PTEN','PTPN11','RAC1','RAF1','RB1','RET',
  'RHOA','RIT1','ROS1','RRAS2','RXRA','SETD2','SF3B1','SMAD3',
  'SMAD4','SMARCA4','SMARCB1','SOS1','SPOP','STAT3','STK11',
  'STK19','TCF7L2','TERT','TGFBR1','TGFBR2','TP53','TP63',
  'TSC1','TSC2', 'U2AF1','VHL','XPO1')

final_colnames = c(
  "TumorId", "NormalId", "Chr1", "Pos1", "Chr2", "Pos2", "SV_Type", 
  "Gene1", "Gene2",  "Transcript1", "Transcript2", "Site1Description", 
  "Site2Description", "Fusion",  "ProbabilityScore", "Confidence", 
  "Comments", "Connection_Type", "SV_LENGTH", "MAPQ",  "PairEndReadSupport", 
  "SplitReadSupport", "BrkptType", "ConsensusSequence",  "TumorReferenceCount", 
  "TumorSplitReferenceCount", "TumorVariantCount",  "TumorSplitVariantCount", 
  "TumorReadCount", "TumorGenotypeQScore", "NormalReferenceCount", 
  "NormalSplitReferenceCount", "NormalVariantCount", "NormalSplitVariantCount", 
  "NormalReadCount", "NormalGenotypeQScore", "Cosmic_Fusion_Counts", 
  "repName-repClass-repFamily:-site1", "repName-repClass-repFamily:-site2", 
  "CC_Chr_Band", "CC_Tumour_Types(Somatic)", "CC_Cancer_Syndrome", "CC_Mutation_Type", 
  "CC_Translocation_Partner", "DGv_Name-DGv_VarType-site1", 
  "DGv_Name-DGv_VarType-site2",  "Significance")

colname_set1 = c(
  "chr1", "pos1", "chr2", "pos2", "gene1", "gene2", 
  "transcript1", "transcript2", "site1", "site2", "fusion", 
  "PairedEndSupport.tumor", "SplitEndSupport.tumor", 
  "PairedEndRef.tumor", "SplitEndRef.tumor", 
  "PairedEndSupport.normal", "SplitEndSupport.normal", 
  "PairedEndRef.normal", "SplitEndRef.normal")

colname_set2 = c(
  "Chr1", "Pos1", "Chr2", "Pos2", "Gene1", "Gene2",
  "Transcript1", "Transcript2", "Site1Description", "Site2Description",
  "Fusion", "TumorVariantCount", "TumorSplitVariantCount",
  "TumorReferenceCount", "TumorSplitReferenceCount", "NormalVariantCount",
  "NormalSplitVariantCount", "NormalReferenceCount", "NormalSplitReferenceCount")


# helper function to parse and retrieve
#  tumor and normal pairedend and split read counts
get_PR.SR_depth = function(samples, data, PRcols=2, SRcols=2) {
  l=list()
  for(sample in samples) {
    for(i in 1:PRcols) {
      for(j in 1:SRcols) {
        l = append(l, list(tstrsplit(
          tstrsplit(data[, get(sample)], ":")[[i]], ",")[[j]]))
      }
    }
  }
  return(l)}


if (!interactive()) {
  parser=ArgumentParser()
  parser$add_argument('-t', '--tab', type='character', help='file name of tab file to be annotated')
  parser$add_argument('-v', '--vcf', type='character', help='file name of vcf file to be annotated')
  parser$add_argument('-o', '--output', type='character', help='file name of output file')
  args=parser$parse_args()
  
  vcf = data.table::fread(args$vcf, colClasses = "character", sep="\t", header=F)
  tab_data = data.table::fread(args$tab, colClasses = "character", sep="\t", header=T)
  
  # Separate vcf metadata
  vcf_meta = vcf[grepl("^#", vcf$V1)]$V1
  vcf_main = vcf[!grepl("^#", vcf$V1)]$V1
  
  # Write header to output file
  write.table(paste(final_colnames, collapse = "\t"),
              args$output, quote = F, row.names = F, 
              col.names = F, append = F)

  # if at least one variant present in vcf
  if (length(vcf_main) > 0) {
    vcf_data = as.data.table(tstrsplit(vcf_main, "\t", fixed = T))
    names(vcf_data) = strsplit(last(vcf_meta), "\t", fixed = T)[[1]]
    
    # Generate genomic coordinates, connection type,
    #  sv length, sv type, and breakpoint type info by
    #  parsting the vcf INFO and FORMAT columns
    vcf_data[, c(
      'chr1', 'pos1', 'chr2', 'pos2', 
      'Connection_Type', 'SV_LENGTH', 'SV_Type', 'BrkptType') := list(
        `#CHROM`, as.character(POS),
        as.character(sub(".*CHR2=([0-9XYM]+).*", "\\1", INFO)),
        as.character(sub(".*END=([1-9]{1}[0-9]+).*", "\\1", INFO)),
        sub(".*CT=([0-9]to[0-9]).*", "\\1", INFO),
        ifelse(grepl("SVLEN", INFO), as.character(
          sub(".*SVLEN=-*([0-9]+).*", "\\1", INFO)), "0"),
        ifelse(sub(".*SVTYPE=([A-Z]+).*", "\\1", INFO) == "BND", "TRA",
               sub(".*SVTYPE=([A-Z]+).*", "\\1", INFO)),
        # set breakpoint type as imprecise if already indicated
        ifelse(grepl('IMPRECISE', INFO), 'IMPRECISE',
        # else, set as precise only if split read support is present
               ifelse(grepl('SR', FORMAT), 'PRECISE', 'IMPRECISE'))),
      by = 1:NROW(vcf_data)]
    
    # Drop the info column
    vcf_data = vcf_data[, !"INFO"]
    
    # Retrieve tumor and normal sample names
    normal_ID = colnames(vcf_data)[9]
    tumor_ID = colnames(vcf_data)[10]
    
    # Parse and extract pairedend and split read info
    vcf_data[, c(
      "PairedEndRef.tumor", "PairedEndSupport.tumor", 
      "SplitEndRef.tumor", "SplitEndSupport.tumor",
      "PairedEndRef.normal", "PairedEndSupport.normal", 
      "SplitEndRef.normal", "SplitEndSupport.normal") := 
        get_PR.SR_depth(c(tumor_ID, normal_ID), vcf_data)]
    
    # Merge vcf table witth the annotation data
    vcf_data = merge(vcf_data, tab_data, 
                     by.x = c("chr1", "pos1", "chr2", "pos2"), 
                     by.y = c("chr1", "pos1", "chr2", "pos2"), 
                     all.x=T, all.y=F)
    
    # Rename columns
    data.table::setnames(vcf_data, colname_set1, colname_set2)
    
    # Add columns for Tumor and Normal samples names, 
    #  pairedend and split read support and total counts,
    #  and significance of the SV based on sig_gene_list list
    vcf_data[, c("TumorId", "NormalId", "PairEndReadSupport", "SplitReadSupport",
                  "TumorReadCount", "NormalReadCount", "Significance") := list(
                    tumor_ID, normal_ID,
                    ifelse(is.na(TumorVariantCount), "0", TumorVariantCount),
                    ifelse(is.na(TumorSplitVariantCount), "0", TumorSplitVariantCount),
                    as.character(
                      sum(as.numeric(TumorVariantCount), 
                          as.numeric(TumorSplitVariantCount),
                          as.numeric(TumorReferenceCount),
                          as.numeric(TumorSplitReferenceCount), na.rm = T)),
                    as.character(
                      sum(as.numeric(NormalVariantCount),
                          as.numeric(NormalSplitVariantCount),
                          as.numeric(NormalReferenceCount),
                          as.numeric(NormalSplitReferenceCount), na.rm = T)),
                    ifelse(Gene1 %in% sig_gene_list | Gene2 %in% sig_gene_list, "KeyGene", "-")), 
              by = 1:NROW(vcf_data)]
    
    # Filter final variants
    #  1. At least one of the gene partner is in ACCESS gene panel
    #  2. At least one of the chromsome partner is an autosome or an allosome
    #  3. Non-translocation SVs should be of 500bp or longer in length
    vcf_data = vcf_data[
      # filter 1
      (vcf_data$Gene1 %in% access_gene_list | vcf_data$Gene2 %in% access_gene_list) &
      # filter 2
      !(vcf_data$Chr1 == 'MT' & vcf_data$Chr2 == 'MT') & 
      # filter 3  
        (SV_Type == 'TRA' | 
           (!SV_Type == "TRA" & as.numeric(SV_LENGTH) >= 500)), ]

    # For BND translocation, with identifical breakpoints, select 3'to5'
    #  connection type to remove duplicates. If connection types are identical
    #  between two variants, select the first. 
    vcf_data[, c("sorted_bkp1", "sorted_bkp2") := list(
	min(paste0(Chr1, ":", Pos1), paste0(Chr2, ":", Pos2)),
	max(paste0(Chr1, ":", Pos1), paste0(Chr2, ":", Pos2))),
    by=1:NROW(vcf_data)]
    
    filter_by_CT = vcf_data[, min(Connection_Type), on = .(Connection_Type), 
	     by = .(SV_Type, sorted_bkp1, sorted_bkp2)]
	     
    vcf_data = unique(vcf_data[filter_by_CT,
	   on = c(SV_Type = "SV_Type", sorted_bkp1 = "sorted_bkp1", 
		  sorted_bkp2 = "sorted_bkp2", Connection_Type = "V1"), 
	   nomatch = 0, mult = "first"])
	     
    
    # Add dummy columns expected to be present in order 
    #  to use the same UI interface as IMPACT
    vcf_data[, c("Confidence", "Comments", "MAPQ", "ConsensusSequence",
                  "TumorGenotypeQScore", "NormalGenotypeQScore",
                  "ProbabilityScore") := rep(list(""), 7)]
    
    # Select for final required columns and write data
    vcf_data = vcf_data[, ..final_colnames]
    write.table(vcf_data, args$output, sep = "\t", row.names = F,
                col.names = F, quote = F, append = T)
  }
}
