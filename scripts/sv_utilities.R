#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
})

# Constants
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
  "DGv_Name-DGv_VarType-site2",  "Significance", "DMPCount", "SplitReadAF")

# Helper functions
resolve_breakpoint = function(x) {
  #' resolve_breakpoint takes a description of a breakpoint of an SV annotated by iAnnotateSV and extracts the region information
  #' @param x (character) iAnnotateSV description of a breakpoint
  #' @return region information of the breakpoint. Example, exon 7, intro 5, 5-UTR

  tryCatch({
  return(ifelse(
    grepl("^Exon", x), tolower(sub(".*(Exon [0-9]+).*", "\\1", x)), ifelse(
      grepl("^Intron", x) & grepl("before exon", x), paste0(
        "intron ", as.character(as.numeric(
          sub(".*before exon ([0-9]+).*", "\\1", x)) - 1)), ifelse(
        grepl("^Intron", x) & grepl("after exon", x), paste0(
        "intron ", sub(".*after exon ([0-9]+).*", "\\1", x)), ifelse(
          grepl("UTR|Promoter|IGR", x), 
          sub("\'", "", sub(
            ".*(5-UTR|3-UTR|5'-UTR|3'-UTR|Promoter|IGR).*", "\\1", x)))))))},
  error = function(e) {
    writeLines(paste("Error when processing breakpoint: ", x, "\n", e)); return("")},
  warning = function(w) {
    writeLines(paste("Warning when processing breakpoint: ", x, "\n", w)); return("")}
  )
}

check_columns = function(sv.dt) {
  #' check and add expected column names in intermediary vcf if necessary
  #' @param sv.dt (data.table) intermediary data table of vcf data
  #' @return vcf data table with expected column names
  #' 
  sv.dt = tryCatch({
    return(sv.dt[, c("EventInfo", "Gene1", "Gene2", 
                     "EventType", "Gene1desc", "Gene2desc")])
  },
  error = function(e) {
    if (grepl("column.*not found", e$message)) {
      sv.dt[, c("EventInfo", "EventType", 
                "Gene1desc", "Gene2desc")] = sv.dt[, c(
                  "Fusion", "SV_Type", 
                  "Site1Description", "Site2Description")]
      return(sv.dt)
    }
    else {
      stop(sprintf('Got this error: %s',e$message))
    }
  })
}

intermediary_columns = function(sv.dt) {
  #' add intermediary columns to vcf data table to describe the genes, breakpoints, and rearrangement genes
  #' @param sv.dt (data.table) intermediary data table of vcf data
  #' @return data table of vcf data with newly added columns

  sv.dt = check_columns(sv.dt)
  sv.dt[, "EVENT" := ifelse(
    grepl("protein fusion", tolower(EventInfo)), "FUSION", ifelse(
      Gene1 == Gene2, "INTRAGENIC", "INTERGENIC")), by=seq_len(NROW(sv.dt))]
  
  sv.dt[, c("GENE1", "GENE2", "SV_GENES") := list(
    ifelse(grepl("TRA", EventType), Gene2, Gene1),
    ifelse(grepl("TRA", EventType), Gene1, Gene2),
    ifelse(EVENT == "FUSION", sub(
      ":", "-", sub(".*[\\(|\\{](.*)[\\)|\\}].*", "\\1", EventInfo)), 
           ifelse(EVENT == "INTRAGENIC", Gene1, 
                  ifelse(grepl("TRA", EventType), 
                         paste0(Gene2, "-", Gene1),
                         paste0(Gene1, "-", Gene2))))),
    by=seq_len(NROW(sv.dt))]
  
  sv.dt[, c("BKP1", "BKP2") := list(
    ifelse(grepl("TRA", EventType), 
           resolve_breakpoint(Gene2desc), resolve_breakpoint(Gene1desc)), 
    ifelse(grepl("TRA", EventType),
           resolve_breakpoint(Gene1desc), resolve_breakpoint(Gene2desc))),
    by=seq_len(NROW(sv.dt))]
  return(sv.dt[,c("EVENT", "GENE1", "GENE2", "SV_GENES", "BKP1", "BKP2")])
}
  
    
#############

make_reference_data = function(
  CohortSVs.f, write.summary=NULL) {
  #' generate summary reference data table of rearrangements and their occurrence
  #' @param CohortSVs.f (character) path to file of reference data
  #' @param write.summary (character) write summary reference data to this file
  #' @return summary reference data table

  if(is.null(raw.cohort.data.f)) {
    stop("Required raw.cohort.data.f argument not defined.")
  }
  CohortSVs = as.data.table(read.csv(
    CohortSVs.f, sep="\t", header=T, colClasses = "character"))
  
  # select minimal essential columns
  CohortSVs = unique(CohortSVs[, c(
    "EventType", "Chr1", "Pos1", "Chr2",     
    "Pos2", "Gene1", "Gene2", "Gene1desc",
    "Gene2desc", "EventInfo")])
  
  CohortSVs = intermediary_columns(CohortSVs)
  CohortSVs = CohortSVs[!(SV_GENES == ""),]
  CohortSVs[, "GENERAL_COUNT" := .N, by=c("SV_GENES")]
  CohortSVs[, "SPECIFIC_COUNT" := .N, by=c("GENE1", "BKP1", "GENE2", "BKP2")]
  
  summary_table = CohortSVs[,c("EVENT", "SV_GENES", "GENE1", 
                                   "BKP1", "GENE2", "BKP2",
                                   "GENERAL_COUNT", "SPECIFIC_COUNT")]
  if(!is.null(write.summary)) {
    write.table(summary_table,
                file=write.summary,
                sep="\t", row.names = FALSE, quote = FALSE)
  }
  return(summary_table)
}

annotate_sv_count = function(
  sv.data.f, write.annotated.f, summary.data.f=NULL, raw.data.f=NULL) {
  #' annotate a file of SVs with occurrence count where applicable
  #' @param sv.data.f (character) path to file of sv data to be annotated
  #' @param write.annotated.f (character) write annotated sv data to this file
  #' @param summary.data.f (character) path to file summary reference data that contains the occurrence count
  #' @param raw.data.f (character) if summary.data.f cannot be provided, a raw reference file can be provided instead. This file will be processed to generate summary reference data with occurrence
  #' @return

  if(is.null(summary.data.f)) {
    !is.null(raw.data.f) || stop(
      "At least one of summary SV data file or a reference file to generate a summary SV data file is required.")
    summary.data = make_reference_data(summary.data.f, write.summary)
  }
  else {
    summary.data = fread(summary.data.f, colClasses = "character", 
                         sep="\t", header = TRUE)
  }

  sv.data = as.data.table(read.csv(sv.data.f, sep="\t",
    header=T, colClasses = "character"))
  cols.to.keep = c(names(sv.data), "DMPCount")
  sv.data = cbind(sv.data, intermediary_columns(sv.data))
  if (nrow(sv.data) > 0) {
    sv.data.intragenic = merge(
      sv.data[(EVENT=="INTRAGENIC")], 
      unique(setnames(summary.data[,!c("GENERAL_COUNT")], 
                      "SPECIFIC_COUNT", "DMPCount")),
      by.x = c("GENE1", "GENE2", "SV_GENES", "BKP1", "BKP2"),
      by.y = c("GENE1", "GENE2", "SV_GENES", "BKP1", "BKP2"),
      all.x = T, no.dups = F)
    
    sv.data.rest = merge(
      sv.data[!(EVENT=="INTRAGENIC")], 
      unique(setnames(summary.data[, c(
        "GENE1", "GENE2", "SV_GENES", "GENERAL_COUNT")],
        "GENERAL_COUNT", "DMPCount")),
      by.x = c("GENE1", "GENE2", "SV_GENES"),
      by.y = c("GENE1", "GENE2", "SV_GENES"),
      all.x = T, no.dups = F)
      
    sv.data = rbind(
      sv.data.intragenic[,..cols.to.keep],
      sv.data.rest[,..cols.to.keep])
  }
  else {
    sv.data[, "DMPCount" :=character()]
    sv.data = sv.data[,..cols.to.keep]
  }

  sv.data[, "SplitReadAF" := round(
    as.numeric(TumorSplitVariantCount)/
      sum(as.numeric(TumorSplitVariantCount),
          as.numeric(TumorSplitReferenceCount)), 5),
  by=seq_len(NROW(sv.data))]
  
  sv.data[is.na(DMPCount)]$DMPCount = "0"
  sv.data[!is.finite(SplitReadAF)]$SplitReadAF = 0
  sv.data[, "SplitReadAF" := as.character(SplitReadAF)]
  
  Chromosomes = c(as.character(1:22), "X", "Y")

  # Reorder rows
  sv.data$Chr1 = factor(sv.data$Chr1, levels=Chromosomes)
  sv.data$Chr2 = factor(sv.data$Chr2, levels=Chromosomes)
  sv.data = sv.data[order(TumorId, Chr1, as.numeric(Pos1),
                          Chr2, as.numeric(Pos2))]
  
  # Rename columns and write table to file
  setnames(sv.data, old = names(sv.data), new = final_colnames)
  write.table(sv.data, write.annotated.f, sep = "\t", row.names = F,
              col.names = T, quote = F, na = "")
}

if (!interactive()) {
  parser = ArgumentParser()
  parser$add_argument('-a', '--annotate-file', type='character', default=NULL, help='file name of structural variants that needs to be annotated with occurrence count')
  parser$add_argument('-o', '--outfile', type='character', default=NULL,
                      help='Name of the output file. If none provided, the input file will be over-written')
  parser$add_argument('-c', '--occurrence-count', type='character', default=NULL,
                      help='Occurrence count data file')
  parser$add_argument('-r', '--raw-cohort-data', type='character', default=NULL,
                      help='raw cohort file that should be used to generate occurrence count data')
  parser$add_argument('-m', '--make-occ-count-data', action="store_true", default = FALSE,
                      help = 'Generate occurrence data. If false, annotate new data')
  
  args = parser$parse_args()
  
  annotate.f = args$annotate_file
  output.f = annotate.f
  if(!is.null(args$outfile)) {
    output.f = args$outfile
  }
  occurrence.count.f = args$occurrence_count
  raw.cohort.data.f = args$raw_cohort_data
  
  if(args$make_occ_count_data) {
    make_reference_data(raw.cohort.data.f, occurrence.count.f)
  }
  else {
   annotate_sv_count(annotate.f, output.f, occurrence.count.f, raw.cohort.data.f)
  }
}
  
  
