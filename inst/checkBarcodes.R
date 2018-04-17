#!/usr/bin/env Rscript
#
#$ -cwd
#$ -S /groups/vbcf-ngs/bin/lib/R-3.2.1/bin/Rscript 

library("forskaller")
library("optparse")
library("plyr")


# bamfile + apikey => create dualbarcodesfile, create splitfile; test 
# dualbarcodesfile + apikey => create splitfile; test
# dualbarcodesfile + splitfile => test
#
#

DUALBC_SCRIPT <- "/groups/vbcf-ngs/bin/funcGen/dualIndexBarcodes.py"
ILLUMINAp5 <- "TCTTTCCCT" 
QIAGENMIRNA <- "GTTCAGAGT"
DEMULTIPLEXER <- "/groups/vbcf-ngs/bin/pipeline/demultiplexer.jar"

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
                            help="complete path to output directory. has to be there already"),
           make_option(c("-n", "--basename"), action="store", dest="basename",
                           help="base name of file to create"),
           make_option(c("-k", "--apikey"), action="store", dest="apikey",
                          help="forskalle 3 API-Key")

)


createDualBarcodesFile <- function(bamfile, dualbarcodesfile){
  res <- system2(DUALBC_SCRIPT, c(bamfile, dualbarcodesfile), stderr=TRUE) 
  if(res != 0){
    error <- paste("error summarising dual barcodes for", bamfile)
    write(error, stderr())
    stop(error)  
  }
}



correctBarcodes <- function(bar){
   toCorrect <- bar$adaptor_tag != "" &  bar$adaptor_secondary_tag == "" & grepl("default:trueseq|illumina", tolower(bar$adaptor_type))
   bar$adaptor_secondary_tag_original <- bar$adaptor_secondary_tag  
   bar$adaptor_secondary_tag[toCorrect] <- ILLUMINAp5 
   bar
}

getAndWriteBarcodes <- function(flowcell, lane, outpath){
   bar <- getBarcodes(flowcell, lane)      
   barc <- correctBarcodes(bar)
   writeBarcodesToFile(barc, outpath) 
   barc
}

check <- function(dualbarcodesfile, splitfile, minLength1, minLength2, prefix){
    cmd <- c("java", " -jar ", DEMULTIPLEXER, " -b ", splitfile, " --len1 ", minLength1, " --len2 ", minLength2, " --mut1 0 --mut2 0 --keepconflicting --check check  ", " -i ", dualbarcodesfile, " -p ", prefix)
    print(cmd)              
    system2(cmd)

}

newRow <- function(bc1, bc2, count, id, bcCount){
   tibble(bc1=bc1, bc2=bc2, count=count, seq=paste(bc1,bc2,sep=""),M1=0,M2=0,outname="",id=id, binnedBarcodes=as.integer(bcCount))
}

reduceBarcodes <- function(tab){
   acc <- tab[1,]
   if(is.na(acc$id)){
     acc <- newRow("", "", acc$count, NA, 1) 
   }
   for(i in 2:nrow(tab)){
      accLast <- acc[nrow(acc),]
      current <- tab[i,]
      if(is.na(accLast$id) && is.na(current$id)){
         sr <- nrow(acc)
         acc$count[sr] <- sum(accLast$count, current$count)
         acc$binnedBarcodes[sr] <- as.integer(acc$binnedBarcodes[sr] + 1)       
         acc$bc1[sr] <- "BINNED" 
         acc$bc2[sr] <- "BINNED"
         acc$seq[sr] <- "BINNED"
      } else {
         acc <- rbind(acc, current)       
      }
   }
   acc 
}



## todo: checkfile path from flowcell/lane

## we return: all to the last identified barcodes max(which(! is.na(checkTable$id))) 
## in extreme cases e.g. wrong barcodes this can be the whole list if the found barcode
## is really at the bottom, so we transform it: 
## 
## only identified barcodes, 
## up to 2x declared barcodes each sequence or minReport(20)
## everything above reportRange (1/10) of last identified barcodes to see strong contaminants, but below 100
## 
## everything unidentified in between only sum of numbers + count of barcodes e.g. 5 barcodes 123456 reads ...
## everything below count of barcodes
##
## we currently search without mismatches
## one could add mismatched identified barcodes with binnedBarcodes 
## 
mangleBarcodeTable <- function(checkTable, barcodes){
    checkTable$binnedBarcodes <- 1:1
    unknown <- checkTable %>% subset(is.na(id))
    known <- checkTable %>% subset(!is.na(id))    

    # summarise by id
    knownS <- known %>% mutate(id2=id) %>% group_by(id) %>% arrange(desc(count)) %>% 
               purrrlyr::by_slice(function(rows){ newRow(rows$bc1[1],rows$bc2[1],count=sum(rows$count), rows$id2[1], nrow(rows)) }, .collate="rows", .labels=FALSE)  
              %>% arrange(desc(count))

    # add summarised rows    
    combTable <- rbind(knownS,unknown) %>% arrange(desc(count))    

    # report at least so many rows
    minReport <- 20
    # report up to reportRange fraction of counts detected with smallest number of barcode e.g. 1/10 of smallest barcode
    reportRange <- 0.1
    # report maximally
    maxReportRange <- 100
    
    minDetected <- max(which(!is.na(combTable$id))) 
    minDetectedCounts <- combTable$count[minDetected]    
    
    rangeMinDetectedCounts <- combTable$count > minDetectedCounts * 0.1
    maxRange <- min(maxReportRange, max(which(rangeMinDetectedCounts)))
    reportMax <- max(minReport, minDetected, maxRange)
    

    report <- combTable[1:reportMax,]
    restBarcodeCount <- nrow(checkTable) - reportMax
    restSum <- sum(checkTable$count[reportMax:nrow(checkTable)])
    lastRow <- newRow("XXXXXXXX", "XXXXXXXX", count=restSum, "rest", restBarcodeCount)
    red <- reduceBarcodes(report)
    redReport <- rbind(red,lastRow)
    
    if(sum(redReport$count) != sum(checkTable$count)){
       stop(paste("counts differ: ",sum(redReport$count), sum(checkTable$count)))
    }
   
    merged <- merge(report, barcodes, by.x="id", by.y="sample_id", all=TRUE)
     
    merged <- merge(checkTable, barcodes, by.x="id", by.y="sample_id", all=TRUE)                 
 
}



doIt <- function(){
   opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = FALSE)
   if(! is.null(opt$bamfile) & !is.null(opt$dualbarcodesfile)){
     createDualBarcodesFile(opt$bamfile, opt$dualbarcodesfile)
   }
   if(! is.null(opt$apikey) & !is.null(opt$splitfile)){
      useKey(opt$apikey)
      bar <- getAndWriteBarcodes(opt$flowcell, opt$lane, opt$splitfile)
      minLength1 <- max(c(6,min(length(bar$adaptor_tag))))
      minLength2 <- min(length(bar$adaptor_secondary_tag))
      r <- check(opt$dualbarcodesfile, opt$splitfile, minLength1, minLength2, checkfile)     
      if(r == 0){
          checkTable <- readr::read_tsv(checkfile)
          found <- checkTable %>% group_by(id) %>% summarise(id=id[1], sum=sum(count)) 
          checkTableTrunc <- max(which(! is.na(checkTable$id)))  
          missing <- merge(bar, found, by.x=c("adaptor_tag", "adaptor_secondary_tag", by.y=c("bc1","bc2"), all.x=TRUE))         
           
      }      
   }

}




