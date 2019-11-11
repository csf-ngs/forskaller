#!/usr/bin/env Rscript
#
#$ -cwd
#$ -S /groups/vbcf-ngs/bin/lib/R-3.2.1/bin/Rscript 

library("optparse")


# bamfile + apikey => create dualbarcodesfile, create splitfile; test 
# dualbarcodesfile + apikey => create splitfile; test
# dualbarcodesfile + splitfile => test
#
#
DUALBC_SCRIPT <- "/groups/vbcf-ngs/bin/funcGen/dualIndexBarcodes.py"

option_list <- list(
           make_option(c("-b", "--bamfile"), action="store", dest="bamfile",
                            help="path to bam file, leave out if starting from dual barcodes file"),
           make_option(c("-d", "--dualbarcodesfile"), action="store", dest="dualbarcodesfile",
                            help="path to file with counts of dual barcodes of bam, will be created when bamfile is used"),
           make_option(c("-s", "--splitfile"), action="store", dest="splitfile", 
                            help="path to file with barcodes to split bam file, will be created if not given"),
           
           make_option(c("-f", "--flowcell"), action="store", dest="flowcell",
                            help="flowcell name"),
           make_option(c("-l", "--lane"), action="store", dest="lane",
                            help="lane"),
            
           make_option(c("-o", "--outputDir"), action="store", dest="outputDir",
                            help="complete path to output directory. will be created"),
           make_option(c("-n", "--basename"), action="store", dest="basename",
                           help="base name of file to create"),
           make_option(c("-k", "--apikey"), action="store", dest="apikey",
                          help="forskalle 3 API-Key")

)


createDualBarcodesFile <- function(bamfile, dualbarcodesfile){
  res <- system2(DUALBC_SCRIPT, c(bamfile, dualbarcodesfile), stdout=TRUE, stderr=TRUE) 
  if(! is.null(attr(res, "status"))){
     error <- paste("error summarising dual barcodes for", bamfile)
     write(error, stderr())
     stop(error)  
  }
}



## todo: checkfile path from flowcell/lane
doIt <- function(){
   opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = FALSE)
   
   dir.create(opt$outputDir, showWarnings = FALSE, recursive = TRUE) 
  
   outPath <- function(fileName){ paste(opt$outputDir, "/", opt$basename, "_", fileName, sep="") }

   if(is.null(opt$splitfile)){
     if(! is.null(opt$apikey)){
       forskaller::useKey(opt$apikey)
       barfile <- outPath("barcodes.tab")
       print("getting barcodes from fsk3")
       bar <- forskaller::getAndWriteBarcodes(opt$flowcell, opt$lane, barfile)
     } else {
       stop("if no splitfile (i.e. barcodes => id mapping) then api is necessary to get it from forskalle3")
     }
   } else {
       bar <- readr::read_tsv(opt$splitfile)
       barfile <- opt$splitfile
   }
   if(is.null(opt$dualbarcodesfile)){
      if(! is.null(opt$bamfile)){
         dualbc <- outPath("all_dual_indices.tab")
         print("extracting indices from bam file")
         createDualBarcodesFile(opt$bamfile, dualbc)
      } else {
         stop("if no dual barcodes file given, must supply bam file to generate it")
      }
   } else {
       dualbc <- opt$dualbarcodesfile
   }
   bar$sample_id <- as.character(bar$sample_id)
   minLength1 <- max(c(6,min(length(bar$adaptor_tag))))
   minLength2 <- min(length(bar$adaptor_secondary_tag))
   checkfile <- outPath("barcode_check.tab")
   r <- forskaller::check(barfile, dualbc, minLength1, minLength2, checkfile)     
   if(r == 0){
       checkTable <- readr::read_tsv(checkfile)
       report <- forskaller::mangleBarcodeTable(checkTable, bar)          
       sumstats <- forskaller::summaryStats(report)      
       write_tsv(report, outPath("barcode_report.tab")) 
       write_tsv(sumstats, outPath("barcode_report_summary.tab"))
       jsonlite::write_json(list(metadata=list(program="barcodecheck", version="0.2", flowcell=opt$flowcell, lane=opt$lane), stats=sumstats, data=report, outPath("barcode_report.json", auto_unbox=TRUE)))
   }      
}

doIt()







