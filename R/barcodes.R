#' @importFrom magrittr %>%


ILLUMINAp5 <- "TCTTTCCCT"
QIAGENMIRNA <- "GTTCAGAGT"
DEMULTIPLEXER <- "/groups/vbcf-ngs/bin/pipeline/demultiplexer.jar"

correctBarcodes <- function(bar){
   toCorrect <- bar$adaptor_tag != "" &  bar$adaptor_secondary_tag == "" & grepl("default:trueseq|illumina", tolower(bar$adaptor_type))
   bar$adaptor_secondary_tag_original <- bar$adaptor_secondary_tag
   bar$adaptor_secondary_tag[toCorrect] <- ILLUMINAp5
   bar
}


#' write barcodes to file
#'
#' @export
getAndWriteBarcodes <- function(flowcell, lane, outpath){
   bar <- getBarcodes(flowcell, lane)
   barc <- correctBarcodes(bar)
   writeBarcodesToFile(barc, outpath)
   barc
}

generateDemultiplexing <- function(barcodes){
  if(nrow(barcodes > 40)){
    mut <- 0
  } else {
    mut <- 1
  }
  b1Min <- min(nchar(barcodes$bc1))
  b1Max <- max(nchar(barcodes$bc1))
  b2Min <- min(nchar(barcodes$bc2))
  b2Max <- max(nchar(barcodes$bc2))
  illumina <- nchar(barcodes$bc1) == 6 & nchar(barcodes$bc2) == 0
  dual <- nchar(barcodes$bc1) > 6 & nchar(barcodes$bc2) > 6
  if(all(illumina)){
    length1 <- b1Min
    length2 <- 0
  } else if(all(dual)){
    length1 <- b1Min
    length2 <- b2Min
  }
  } else if(any(illumina) & any(dual)){
    length1 <- b1Min
    length2 <- 8
  } else if(! any(illumina) & all(dual)){
    length1 <- b1Min
    length2 <- b2Min
  } 

  c( "--mut1", mut, "--mut2", mut,  "--len1", minLength1, "--len2", minLength2 ) 
  

}



#' check barcodes file
#'
#' @export
check <- function(dualbarcodesfile, splitfile, prefix, otherOptions){
    cmd <- c("java", "-jar", DEMULTIPLEXER, "-b", dualbarcodesfile, otherOptions, "--keepconflicting", "--check", "check", "-i", splitfile, "--outpathcheck", prefix)
    print(paste(cmd, collapse=" "))
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

getBarcodeOfLength <- function(char, len){
   paste(rep(char, len), collapse="")
}

#' make summary stats
#'
#' @export
summaryStats <- function(report){
   report %>% dplyr::group_by(!is.na(id) & id != "NNNN" & id != "BINNED") %>% dplyr::group_by(multi_id) %>% dplyr::summarize(count=sum(count)) %>% rename_at(1, ~"known") %>% mutate(ratio=count/sum(count))
}

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
#' mangle barcode table 
#'
#' @export
mangleBarcodeTable <- function(checkTable, barcodes){
    checkTable$binnedBarcodes <- 1:1
    unknown <- checkTable %>% subset(is.na(id))
    known <- checkTable %>% subset(!is.na(id))

    # summarise by id
    knownS <- known %>% mutate(id2=id) %>% group_by(id) %>% arrange(desc(count)) %>%
               purrrlyr::by_slice(function(rows){ newRow(rows$bc1[1],rows$bc2[1],count=sum(rows$count), rows$id2[1], nrow(rows)) }, .collate="rows", .labels=FALSE) %>%
               arrange(desc(count))

    # add summarised rows
    combTable <- rbind(knownS,unknown) %>% arrange(desc(count))

    # report at least so many rows
    minReport <- 20
    # report up to reportRange fraction of counts detected with smallest number of barcode e.g. 1/10 of smallest barcode
    reportRange <- 0.1
    # report maximally
    maxReportRange <- 100

    bcLength1 = nchar(checkTable[1,]$bc1)
    bcLength2 = nchar(checkTable[1,]$bc2)
    Nbc1 = getBarcodeOfLength("N", bcLength1)
    Nbc2 = getBarcodeOfLength("N", bcLength2)
    Xbc1 = getBarcodeOfLength("X", bcLength1)
    Xbc2 = getBarcodeOfLength("X", bcLength2)

    minDetected <- max(which(!is.na(combTable$id)))
    minDetectedCounts <- combTable$count[minDetected]

    NNNrow <- which(combTable$bc1 == Nbc1 & combTable$bc2 == Nbc2)
    if(length(NNNrow) > 0){
      combTable$id[NNNrow] <- "NNNNN"
    }

    rangeMinDetectedCounts <- combTable$count > minDetectedCounts * 0.1
    maxRange <- min(maxReportRange, max(which(rangeMinDetectedCounts)))
    reportMax <- max(minReport, minDetected, maxRange)
    reportTop <- combTable[1:reportMax,]
    reportBottom <- combTable[-(1:reportMax),]
    reportBottomRed <- reduceBarcodes(reportBottom)

    report <- rbind(reportTop, reportBottomRed)

    if(sum(report$count) != sum(checkTable$count)){
       stop(paste("counts differ: ",sum(redReport$count), sum(checkTable$count)))
    }

    merged <- full_join(report, barcodes, by=c("id"="sample_id"))
    merged
}



