suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
})

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
          grepl("5-UTR|3-UTR|Promoter|IGR", x), 
          sub(".*(5-UTR|3-UTR|Promoter|IGR).*", "\\1", x))))))},
  error = function(e) {
    writeLines(paste("Error when processing breakpoint: ", x, "\n", e)); return("")},
  warning = function(w) {
    writeLines(paste("Warning when processing breakpoint: ", x, "\n", w)); return("")}
  )
}

make_reference_data = function(
  raw.data, write.summary=NULL) {
  signedout_svs = as.data.table(read.csv(
    "/Volumes/dmp-hot/zehira/All_Results/All_so_IMPACT_Solid_SVs.txt", sep="\t",
    header=T, colClasses = "character"))
  
  # select minimal essential columns
  signedout_svs = signedout_svs[, c(
    "EventType", "Chr1", "Pos1", "Chr2",     
    "Pos2", "Gene1", "Gene2", "Gene1desc",
    "Gene2desc", "EventInfo")]
  
  signedout_svs = unique(signedout_svs)
  signedout_svs = signedout_svs[, !c("DMP_SAMPLE_ID")]
  signedout_svs = signedout_svs_bkp
  
  signedout_svs[, "EVENT" := ifelse(
    grepl("protein fusion", tolower(EventInfo)), "FUSION", ifelse(
      Gene1 == Gene2, "INTRAGENIC", "INTERGENIC")), by=1:NROW(signedout_svs)]
  
  signedout_svs[, c("GENE1", "GENE2", "SV_GENES") := list(
    ifelse(grepl("TRA", EventType), Gene2, Gene1),
    ifelse(grepl("TRA", EventType), Gene1, Gene2),
    ifelse(Event == "FUSION", sub(".*\\((.*)\\).*", "\\1", EventInfo), 
           ifelse(Event == "INTRAGENIC", Gene1, 
                  ifelse(grepl("TRA", EventType), 
                         paste0(Gene2, "-", Gene1),
                         paste0(Gene1, "-", Gene2))))),
    by=1:NROW(signedout_svs)]
  
  signedout_svs[, c("BKP1", "BKP2") := list(
    ifelse(grepl("TRA", EventType), 
           resolve_breakpoint(Gene2desc), resolve_breakpoint(Gene1desc)), 
    ifelse(grepl("TRA", EventType),
           resolve_breakpoint(Gene1desc), resolve_breakpoint(Gene2desc))),
           by=1:NROW(signedout_svs)]
  
  signedout_svs = signedout_svs[!(SV_GENES == ""),]
  signedout_svs[, "GENERAL_COUNT" := .N, by=c("SV_GENES")]
  signedout_svs[, "SPECIFIC_COUNT" := .N, by=c("GENE1", "BKP1", "GENE2", "BKP2")]
  
  summary_table = signedout_svs[,c("EVENT", "SV_GENES", "GENE1", 
                                   "BKP1", "GENE2", "BKP2",
                                   "GENERAL_COUNT", "SPECIFIC_COUNT")]
  if(!is.null(write.summary)) {
    write.table(summary_table,
                file=write.summary,
                sep="\t", row.names = FALSE, quote = FALSE)
  }
  return(summary_table)
}


