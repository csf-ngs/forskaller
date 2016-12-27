require("plyr")
require("RCurl")
require("rjson")
require("Hmisc") #latexTranslate
require("httr")
require("stringr")

##
## http://ngs.vbcf.ac.at/forskalle/apidoc

#' create credentials
#' returns list username=username password=password
#'
#' @param username your user name
#' @param password your password
#' @export
createCredentials <- function(username, password){
  list(username=username, password=password)
}


#' starts a forskalle session
#' the api is a little stateful
#' you should call endSession at the end
#' returns a forskalle session that subsequent interactions with forskalle need
#'
#' its: session <- startSession(createCredentials(username, password))
#' 
#' @param credentials list created with createCredentials
#' @export
startSession <- function(credentials){
  loginurl <- "http://ngs.vbcf.ac.at/forskalle/api/login"
  curl <- getCurlHandle()
  curlSetOpt(cookiejar="", followlocation = TRUE, curl=curl)
  tryCatch(
    loginResult <- postForm(loginurl, .params = credentials, curl=curl, .checkParams=FALSE ),
    error = function(e) { cat(paste("problem creating session with username: ", credentials$username, "\n", e), file=stderr())} 
  )
  curl
}


#' ends a forskalle session
#'
#' @param session
#' @export
endSession <- function(session){
  logouturl <- "http://ngs.vbcf.ac.at/forskalle/logout" #GET
  getURLContent(logouturl, curl=session)
  rm(session)
}

## replace null with NA for as.data.frame
nullToNA <- function(li){
  Map(function(v){ if(is.null(v)){ NA }else{ v }}, li)
}

### its a list of individual measurements. The form attribute gives the type ###
## mostly its just one value. changesets is a list => remove, but keep id
## form:
## "Preparation" , "Size Analysis", "qPCR" , "RNA Quantification"
measurementToDF <- function(meas){
  ci <- which(names(meas) == "changesets")  
  batchId <- NA
  if(length(meas$changesets) > 0){
     batchId <- meas$changesets[[1]]$id
  }
  mf <- meas[-ci]
  mf$batchId <- batchId
  mfn <- nullToNA(mf)
  list(type=mfn$form, data=as.data.frame(mfn, stringsAsFactors=FALSE))
}

removeFromList <- function(lis, re){
  li <- which(names(lis) == re)
  lis[-li]
}


#' gets sample information from forskalle
#'
#' @param sampleId the sample id 
#' @param session the forskalle session
#' @param should measurement columns be simplified  
#'
#' returns a list of
#' sample = data frame with sample info (1 row)
#' measurments = a list of data frames with measurements 
#'
#' @export
getSample <- function(sampleId, session, simplify=FALSE){
  s <- NULL
  tryCatch(
   s <- getURLContent(paste("http://ngs.vbcf.ac.at/forskalle/api/samples/", sampleId, sep=""), curl=session), ## its a string,
   error=function(e){ cat(paste("error retrieving sample info: ", sampleId, "\n", e), file=stderr()) }
  )
  if(is.null(s)){
    return(s)
  }
  sj <- fromJSON(s) ## its a nested list
  sjl <- removeFromList(sj, "requests")
  #logicalColumns <- c("shearing", "add_primer", "own_risk", "fragmented", "stranded")
  logicalColumns <- c("shearing", "own_risk", "fragmented", "stranded")
  logicalColumnsIndex <- names(sjl) %in% logicalColumns 
  sjl[logicalColumnsIndex] <- as.logical(sjl[logicalColumnsIndex])
  # must replace null with NA
  sjln <- nullToNA(sjl)
  
  sampleDF <- as.data.frame(sjln,stringsAsFactors=FALSE)
  sampleDF <- rename(sampleDF, c("preparation_type"="prep"))
  sampleDF$id <- as.integer(sampleDF$id) # is numeric ?? 
  tryCatch(
     m <- getURLContent(paste("http://ngs.vbcf.ac.at/forskalle/api/measurements/sample/", sampleId, sep=""), curl=session), ## its a string    
     error=function(e){ cat(paste("error retrieving measurement info: ", sampleId, "\n", e), file=stderr()) }
  )
  mj <- fromJSON(m)
  mea <- lapply(mj, measurementToDF)
  if(simplify){
     sm <- lapply(mea, simplifyMeasurement)     
     sa <- simplifySample(sampleDF)
     #sample table is very long, split it somehow in lab/annotation
     #"id"               "preparation_type" "cutout_size"      "shearing"         "fragmented"       "stranded"         "own_risk"         "add_primer"       "exptype"          "organism"         "genotype"         "celltype"         "antibody"         "descr"
     sampleLab <- subset(sa, select=c(id, prep, cutout_size, shearing, fragmented, stranded, own_risk))#, add_primer))
     sampleAnnot <- subset(sa, select=c(id, exptype, organism, genotype, celltype, antibody, descr))
     list(sample=sa, sampleLab=sampleLab, sampleAnnot=sampleAnnot, measurements=sm)
  } else {
     list(sample=sampleDF, measurements=mea) 
  }

}

#' gets samples as list of data frames
#' 
#' @param sampleIds vector of sample ids
#' @param session the forskalle session
#'  
#' returns a list of data frames
#' samples
#' for each measurement type a data frame
#' all data frames are simplified (as in getSample)
#'
#' @export
getSamples <- function(sampleIds, session){
   samples <- NULL
   tryCatch(
      samples <- lapply(sampleIds, getSample, session, TRUE),
      error=function(e){ cat(paste("error retrieving samples info: ", sampleIds, "\n", e), file=stderr()) }             
   )
   isNullS <- sapply(samples, function(s){ is.null(s) })
   if(any(isNullS)){
      return(NULL) 
   }
   samplesDF <- do.call("rbind", lapply(samples, function(s){ s$sample }))
   samplesADF <-  do.call("rbind", lapply(samples, function(s){ s$sampleAnnot }))
   samplesLDF <-  do.call("rbind", lapply(samples, function(s){ s$sampleLab }))

   measurementTypesM <- sapply(samples, function(s){  
      sapply(s$measurements, function(m){ m$type })
   })
   measurementTypes <- unique(unlist(measurementTypesM))

   measurements <- lapply(measurementTypes, function(type){
        meas <- Map(function(s){ s$measurements }, samples)
        ms <- Filter(function(m){ m$type == type }, unlist(meas, recursive=FALSE))
        dat <- do.call("rbind", Map(function(m){ m$data }, ms))
        dato <- dat[order(dat$sampleId),] 
        dato
   })
   names(measurements) <- measurementTypes    
   sampleso <- samplesDF[order(samplesDF$id),]
   samplesao <- samplesADF[order(samplesADF$id),]
   sampleslo <- samplesLDF[order(samplesLDF$id),]

   list(samples=sampleso, samplesAnnot=samplesao, samplesLab=sampleslo, measurements=measurements)
}

#' get one multiplex info
#' 
#' 
#' data.frame(multiId, sampleId, tag, ratio)
#'
#' @param multiId the id  multiplex ids start with M TODO: could check this? 
#'                        if no M then some other info comes up
#' @param session the session
#'
#' @export
getMultiplex <- function(multiId, session){
   tryCatch(
     multi <- getURLContent(paste("http://ngs.vbcf.ac.at/forskalle/api/multiplexes/", multiId, sep=""), curl=session), ## its a string,
     error=function(e){ cat(paste("error retrieving multiplex info: ", multiId, "\n", e), file=stderr()) }
   )	
   mj <- fromJSON(multi)
   mjs <- mj$samples
   sb <- do.call("rbind", lapply(mjs, function(s){ data.frame(sampleId=s$sample$id, tag=s$sample$tag, ratio=s$ratio, stringsAsFactors=FALSE)}))
   sb$multiId <- multiId
   subset(sb,select=c(multiId,sampleId,tag,ratio)) 
}

#' compare 2 multiplexes by barcode
#' 
#'
#' @export
compareMultiplex <- function(m1,m2){
	m <- merge(m1, m2, by.x="tag", by.y="tag", all.x=TRUE, all.y=TRUE)
	m$unique <- is.na(m$multiId.x) | is.na(m$multiId.y)	
    mo <- m[order(m$sampleId.x, m$sampleId.y),]
    mo
}


#' selects specific columns from measurement
#' @param measurement (list type=type, data=data)
#'
#' @export
simplifyMeasurement <- function(measurement){
    simple <- switch(measurement$type,
        'Preparation'=simplifyPreparation(measurement$data),
        'Size Analysis'=simplifySizeAnalysis(measurement$data),
        'qPCR'=simplifyQPCR(measurement$data),
        'RNA Quantification'=simplifyRNAQuantification(measurement$data),
        'Quantification'=simplifyQuantification(measurement$data),
        'cDNA Synthesis'=simplifyCDNASynthesis(measurement$data), 
        stop(paste("unknown measurement ", measurement$type))
    )
    sr <- rename(simple, c("obj_id"="sampleId", "multi_id"="multiId"))
    list(type=measurement$type, data=sr)
}

#' subsets but is lenient about missing columns
#' missing columns are set to NA
#' @param df the data frame to subset by columns
#' @param cols string vector of columns to subset in desired order
#'
subsetF <- function(df, cols){
   coli <- match(cols, colnames(df))
   colp <- coli[!is.na(coli)]
   missingColumns <- which(is.na(coli))
   allColumns <- 1:length(cols)
   dfp <- df[,colp]
   alldf <- NULL
   nacol <- rep(NA, nrow(df))
   dfi <- 0
   for(ci in allColumns){
     col <- if(ci %in% missingColumns){ nacol }else{ dfi = dfi + 1; dfp[,dfi] }
     alldf <- cbind(alldf, col)
   }
   colnames(alldf) <- cols
   data.frame(alldf)
}


#' simplifies preparation data.frame
#'
#'   severity    id       type        form         user obj_id cycles udgase text notified cutout_size                date flag  change_user         change_date multi_id obj_type resolved           kit method batchId
#'      0 11344 Data Entry Preparation Carmen Czepe  16864     12      1   NA        0     200-800 2014-02-24 14:00:36   Ok Carmen Czepe 2014-02-26 16:28:39     1182   sample        0 NEB ultra RNA  beads    1391
#'
simplifyPreparation <- function(preparation){
 subsetF(preparation, c("obj_id", "batchId", "multi_id", "cycles", "udgase", "cutout_size", "flag", "user", "kit", "method", "date"))  
}


#' simplifies sample data.frame
#'
#'   tissue_type preparation_type cutout_size obj_id    tag add_primer    id celltype own_risk status   received obj_type shearing         barcode      ready fragment_size       scientist fragmented    group secondary_tag             descr stranded
#'                       RNA           -  16864 ATCACG      FALSE 16864     MEFs    FALSE  Ready 2014-01-14   sample    FALSE Default:TrueSeq 2014-03-03             - Jennifer Jurkin       TRUE Martinez            NA RTCB -OHT d4 E189        0
#'    primer tagno                              comments exptype antibody genotype     organism secondary_tagno userprep
#' 1 Standard     1 libary preparation beginning from RNA RNA-Seq          wildtype Mus musculus              NA 
simplifySample <- function(sample){
  truncateTo <- function(string){
    if(is.na(string) || is.null(string)){
       return("NA")
    }
    stringc <- as.character(string)
    tru <- if(nchar(stringc) > 25){ 
      s <- substr(stringc, 1, 25)
      paste(s, "...", sep="")
    }
    else{
      stringc
    }
    latexTranslate(tru)
  }
  #subs <- subset(sample, select=c(id, tag,  prep, cutout_size, shearing, fragmented, stranded, own_risk, add_primer, exptype, organism, genotype, celltype, antibody, descr))
  subs <- subset(sample, select=c(id, tag,  prep, cutout_size, shearing, fragmented, stranded, own_risk, exptype, organism, genotype, celltype, antibody, scientist, group, descr, preparation_kit))
  subs$description <- subs$descr
  within(subs, { genotype=truncateTo(genotype); celltype=truncateTo(celltype); antibody=truncateTo(antibody); descr=truncateTo(descr)})
}


#' simplify RNA Quantification data.frame
#'
#'   resolved conc rin multi_id obj_type                date text notified         change_date  change_user flag       type               form severity    id obj_id         user batchId
#       0  618 8.6     1182   sample 2014-02-20 15:23:32   NA        0 2014-02-21 15:14:57 Carmen Czepe   Ok Data Entry RNA Quantification        0 11278  16864 Carmen Czepe    1376
#'
simplifyRNAQuantification <- function(quantification){
  subsetF(quantification, c("obj_id", "batchId", "multi_id", "conc", "rin", "flag", "user", "date"))
}

#' simplify size analysis data.frame
#' 
#' id severity dilution          form       type        user obj_id notified text                date flag change_user         change_date multi_id obj_type resolved molarity size kit conc method batchId
#' 1 11436        0    -0.51 Size Analysis Data Entry Laura Bayer  16864        0   NA 2014-02-27 14:14:27   Ok Laura Bayer 2014-02-27 15:22:28     1182   sample        0     4.15  270  FA 0.74     HS    1429
simplifySizeAnalysis <- function(sizeanalysis){
  subsetF(sizeanalysis, c("obj_id", "batchId", "multi_id", "dilution", "size", "conc", "flag", "kit", "method", "user", "date")) 
}


#' simplify qpcr data.frame
#'
#'   resolved size efficiency  kit conc multi_id obj_type X2nM_control corrected_conc notified machine text                date flag change_user         change_date severity    id       type form obj_id     user R2 batchId
#'         0  270       96.9 Kapa  2.9     1182   sample         2.24           4.85        0    ours      2014-03-03 11:44:00   Ok    Ru Huang 2014-03-03 11:47:20        0 11478 Data Entry qPCR  16864 Ru Huang  1    1437 
simplifyQPCR <- function(qpcr){
  subsetF(qpcr, c("obj_id", "batchId", "multi_id", "size", "efficiency", "conc", "corrected_conc", "kit", "flag", "user", "date"))
}

#' simplify Quantification (Chip-Seq)
#'
#'
#'      user obj_type volume multi_id notified         change_date bioanalyzer_result obj_id           form resolved change_user  conc       type                date    id total text flag severity bioanalyzer_done batchId
#'     Ru Huang   sample      9     1250        0 2014-03-04 12:52:36                 Ok  17853 Quantification        0    Ru Huang 10010 Data Entry 2014-03-04 11:10:39 11560 90.09   NA   Ok        0                1    1445
simplifyQuantification <- function(quant){
  subsetF(quant, c("obj_id", "batchId", "multi_id", "conc", "flag", "user", "date"))
}

#' simplify cDNA synthesis (RNA-Seq)
#' 
#'
#'
simplifyCDNASynthesis <- function(cdna){
  subsetF(cdna, c("obj_id", "batchId", "multi_id", "flag", "kit", "user", "ercc", "date"))
}

#' get today formatted for forskalle
#'
today <- function(){
  format(Sys.time(), "%Y-%m-%d")
}

#' create A sample
#'
#'
#' @export
#createSample <- function(description, comments, group, scientist, celltype="", genotype="", exptype="ChIP-Seq", organism="", tissue_type="", antibody="", status="Aborted", received=today(), fragmented=1, userprep=1, preparation_kit=NULL, preparation_type="none", own_risk=0, barcode="", stranded=0,  shearing=0, secondary_tag=NULL, add_primer=0, cutout_size="100-700", tagno=NULL, secondary_tagno=NULL,  ready=NULL, primer="Standard", cutout_size_min=100, cutout_size_max=700, fragment_size=""){
createSample <- function(description, comments, group, scientist, celltype="", genotype="", exptype="ChIP-Seq", organism="", tissue_type="", antibody="", status="Aborted", received=today(), fragmented=1, userprep=1, preparation_kit=NULL, preparation_type="none", own_risk=0, barcode="", stranded=0,  shearing=0, secondary_tag=NULL, cutout_size="100-700", tagno=NULL, secondary_tagno=NULL,  ready=NULL, primer="Standard", cutout_size_min=100, cutout_size_max=700, fragment_size=""){
   li <- list(
      description=description,
      comments=comments,
      group=group,
      scientist=scientist,
      celltype=celltype,
      genotype=genotype,
      exptype=exptype,
      antibody=antibody,
      organism=organism,
      tissue_type=tissue_type,
      status=status,
      received=received,
      fragmented=fragmented,
      userprep=userprep,
      preparation_kit=preparation_kit,
      preparation_type=preparation_type,
      own_risk=own_risk,
      barcode=barcode,
      secondary_tag=secondary_tag,
    #  add_primer=add_primer,
      cutout_size=cutout_size,
      tagno=tagno,
      secondary_tagno=secondary_tagno,
      ready=ready,
      primer=primer,
      cutout_size_min=cutout_size_min,
      cutout_size_max=cutout_size_max,
      fragment_size=fragment_size,
      requests={},
      isNew=TRUE
   )
   li

#{"secondary_tag":null,"status":"Aborted","exptype":"ChIP-Seq","obj_type":"sample","id":26038,"shearing":0,"genotype":"wildtype","descr":"SRR955859 Pol II ChIP resting B cells (Kouzine et al) rep1","group":"Pavri","requests":[],"add_primer":0,"cutout_size":"100-700","tagno":null,"celltype":"B cells","scientist":"Rushad Pavri","organism":"Mus musculus","fragmented":1,"userprep":"1","obj_id":26038,"preparation_kit":null,"preparation_type":"none","own_risk":0,"barcode":"Nextera","stranded":0,"received":"2015-02-23","ready":null,"comments":"SRR955859","tissue_type":"","primer":"Standard","cutout_size_min":"100","antibody":"Pol II","secondary_tagno":null,"cutout_size_max":"700","fragment_size":"","tag":"random"}

}

#' add sample to forskalle
#' @param sample created with createSample
#'
#' @export
addSample <- function(sample, session){
    sampleJson <- toJSON(sample)
    httpheader <- c(Accept="application/json; charset=UTF-8", "Content-Type"="application/json")
    tryCatch(
       addResult <- postForm("http://ngs.vbcf.ac.at/lammskalle/api/samples", curl=session, .opts=list(postfields=sampleJson))
 ,
       error = function(e) { cat(paste("problem creating session with username: ", credentials$username, "\n", e), file=stderr())}
    )
}
 
#' turns one checkresults to a data frame
#'
#' @export
checkResultsToDF <- function(checkResultsIndex, checkResultsList, sampleId){
   checkResults <- checkResultsList[[checkResultsIndex]]
   fp<- checkResults[[1]]
   read2q30 <- NA
# out of bounds but not necessary now
   if(length(checkResults) > 1){
      sr <- checkResults[[2]]
      read2 <- sr$read
      if(read2 == "1"){
         fp <- sr
         read2q30 <- checkResults[[1]]$countq30 
      }else{
         read2q30 <- sr$countq30
      }
   }
   data.frame(sampleId=sampleId, basecallsNr=LETTERS[checkResultsIndex], flowcell=fp$flowcell, lane=fp$lane, basecalls=fp$basecalls, result=fp$result, total=fp$total, q30.1=fp$countq30, q30.2=read2q30)
}

#' get run info from forskalle
#' 
#' @export
runsForSample <- function(sampleId, session){
  s <- NULL
  tryCatch(
   s <- getURLContent(paste("http://ngs.vbcf.ac.at/forskalle/api/runs/sample/", sampleId, sep=""), curl=session), ## its a string,
   error=function(e){ cat(paste("error retrieving run info for sample: ", sampleId, "\n", e), file=stderr()) }
  )
  if(is.null(s)){
    return(s)
  }
  sj <- fromJSON(s) ## its    
  if(length(sj) > 0){
  checkResults <- sj[[1]]$check_results
  allChecks <- do.call("rbind", lapply(seq_along(checkResults), checkResultsToDF, checkResults, sampleId))
  allChecks
 }else{
     data.frame(sampleId=sampleId, basecallsNr=NA, flowcell=NA, lane=NA, basecalls=NA, result=NA, total=NA, q30.1=NA, q30.2=NA)
  }
}



#' get samples for user
#'
#' @export
samplesForGroup <- function(groupName, session){
  s <- NULL
  tryCatch(
   s <- getURLContent(paste("http://ngs.vbcf.ac.at/forskalle/api/samples?group=", groupName, "&since=10-1", sep=""), curl=session), ## its a string,
   error=function(e){ cat(paste("error retrieving samples for group: ", groupName, "\n", e), file=stderr()) }
  )
  if(is.null(s)){
    return(s)
  }
  sj <- fromJSON(s) ## its
  sapply(sj, function(s){ s$id })
}

multiplexToDf <- function(multiplex){
   samples <- multiplex$samples
   mdf <- do.call("rbind", lapply(samples, function(s){
     data.frame(id=s$sample$obj_id, barcode=s$sample$tag, ratio=as.numeric(s$ratio), stringsAsFactors=FALSE) 
   }))
   mdf
}



#' creates long data.frame from multiplex
#'
#' G/T = green laser
#' A/C = red laser
#'
multiplexDfToLongWithBarcodes <- function(multiplexDF){
   ddply(multiplexDF, .(id), function(s){
      base <- strsplit(s$barcode, "")[[1]]
      position <- seq_along(base)
      laser <- ifelse(base == "A" | base == "C", "red", "green")
      data.frame(id=s$id, base=base, position=position, laser=laser, ratio=s$ratio)
   })
}

#' get multiplex 
#'
#'
multiplexById <- function(multiplex, session){
  s <- NULL
  tryCatch(
   s <- getURLContent(paste("http://ngs.vbcf.ac.at/forskalle/api/multiplexes/", multiplex, sep=""), curl=session), ## its a string,
   error=function(e){ cat(paste("error retrieving multiplex: ", multiplex, "\n", e), file=stderr()) }
  )
  if(is.null(s)){
    return(s)
  }
  sj <- fromJSON(s) ## its
  multiplexToDf(sj)
}


#' make barcode chart for multiplex
#' 
#' 
#' 
barcodeChartForLane <- function(multiplexLong){
  ggplot(multiplexLong, aes(x=position, y=as.factor(id), fill=base)) + geom_tile() + ylab("sample id") + scale_fill_manual(values=c("A"="yellow1","C"="yellow2","G"="royalblue1","T"="royalblue2"))  
}

columnBaseCounts <- function(base, bases, counts){
  sum(counts[base == bases])
}

# the X is in the database as XXXXXXX for not in top 100 barcode combinations
baseCounts <- function(bases, counts){
  nucs <- c("A","C","G","T","N","X")
  nc <- lapply(nucs, columnBaseCounts, bases, counts)
  names(nc) <- nucs
  as.data.frame(nc)
}

allBaseCounts <- function(seqsDF, counts){
  bcs <- apply(seqsDF, 2, baseCounts, counts)
  bdf <- do.call("rbind", bcs) 
  xSum <- sum(bdf$X)
  if(xSum == 0){
    bdf[,-ncol(bdf)]    
  } else {
    bdf
  }
}

barcodeCountsToRatio <- function(path){
   tab <- read.table(path, header=TRUE, sep="\t", as.is=TRUE)
   seqs <- strsplit(tab[,1], "")
   seqsDF <- do.call("rbind", seqs)
   counts <- tab[,2]
   basedf <- allBaseCounts(seqsDF, counts)
   rs <- rowSums(basedf)
   freqs <- basedf/rs
   fdf <- data.frame(position=1:nrow(freqs), freqs) 
   fdf
}


#' create from barcode a sequencelogo
#' @param path ???
#'
#'
#' @export
sequenceFreqs <- function(path){
   fdf <- barcodeCountsToRatio(path) 
   fdfl <- melt(fdf, id.vars="position")       
   ggplot(fdfl, aes(x=factor(position),y=value, fill=variable)) + geom_bar(position="stack", stat="identity") + guides(fill=guide_legend("base")) + ylab("ratio") + xlab("index cycle") 
}


#' get todays day for to for getRuns
#' 
#' 
#' @export
getToday <- function(){
   format(Sys.time(), "%Y-%m-%d")
}


#' removes multiple elements from a list by names
#' @param li the list 
#' @param vector of names to remove from list 
#'
#' @export
removeElementsFromList <- function(li, multinames){
    indices <- which(names(li) %in% multinames)
    lic <- li[-indices]
    lic <- lic[order(names(lic))]
}


#' convert flowcell to table row + list
#' @param flowcell
#'
#'
#' @export
flowcellToTable <- function(flowcell){
   toRemove <- c("lanes", "problems")   
   values <- removeElementsFromList(flowcell, toRemove) 
   values 
}

#' gets tag length even if null
#'
#' @param tag
#'
#' @export
getTagLength <- function(tag){
  if(is.null(tag)){ 0 }else{ nchar(tag) }
}

#' get tags length for sample 
#' @param sample
#' 
#' @export
getTagsForSample <- function(sample){
   ss <- sample$request_sample$sample
   barcode <- ss$barcode
   ltag1 <- getTagLength(ss$tag)
   ltag2 <- getTagLength(ss$secondary_tag)
   spike <- if(is.null(sample$is_spikein)){ FALSE }else{ sample$is_spikein == 1 }
   scientist <- ss$scientist
   description <- ss$descr
   data.frame(barcode=barcode, ltag1=ltag1, ltag2=ltag2, spikeIn=spike, scientist=scientist, description=description)
}


#' iterate through samples on lane and get length of tags
#' lane$samples[[1]]$request_sample$sample
#' @param samples
#'
#' @export
getTagsForSamples <- function(samples){
   do.call("rbind", lapply(samples, getTagsForSample))
}

simplifyTagsForLane <- function(tags){
   ltag1 <- max(tags$ltag1)
   ltag2 <- max(tags$ltag2)
   hasSpikeIn <- any(tags$spikeIn)
   scientist <- unique(tags$scientist)
   scientist <- if(length(scientist) > 1){ "multi" }else{ scientist }
   description <- unique(tags$description)
   description <- if(length(description) > 1){ "multi" }else{ description }
   data.frame(ltag1=ltag1, ltag2=ltag2, hasSpikeIn=hasSpikeIn, sampleCount=nrow(tags), scientist=scientist, description=description)
}

#' get simplified tags for lane
#' @param lane
#'
#' @export
getTagsForLane <- function(lane){
   laneNr <- lane$num
   tags <- getTagsForSamples(lane$samples)
   simp <- simplifyTagsForLane(tags)
   data.frame(lane=laneNr, simp)
}

#' get flowcell by id 
#' @param id
#' @param session
#'
#' @export
getFlowcellById <- function(id, session){
  s <- NULL
  query <- paste("http://ngs.vbcf.ac.at/forskalle/api/flowcells/", id, sep="")
  tryCatch(
   s <- getURLContent(query, curl=session), ## its a string,
   error=function(e){ cat(paste("error retrieving  flowcell ", id, "\n", e), file=stderr()) }
  )
  if(is.null(s)){
    return(s)
  }
  sj <- fromJSON(s) ## its a nested list
  laneTags <- do.call("rbind", lapply(sj$lanes, getTagsForLane))
  laneTags$flowcell <- id
  laneTags
}


#' get flowcells  ?from=2014-09-02&to=2015-12-02
#' @param from
#' @param to
#' from <- "2014-09-02"
#'
#' @export
getFlowcells <- function(session, from="2014-01-02", to=getToday()){
  s <- NULL
  query <- paste("http://ngs.vbcf.ac.at/forskalle/api/flowcells?from=",from, "&to=", to, sep="")
  tryCatch(
   s <- getURLContent(query, curl=session), ## its a string,
   error=function(e){ cat(paste("error retrieving flowcells info ", from, "-", to, "\n", e), file=stderr()) }
  )
  if(is.null(s)){
    return(s)
  }
  sj <- fromJSON(s) ## its a nested list
  sjc <- lapply(sj, removeElementsFromList, c("problems", "lanes", "planned_start", "comments"))    
  rbf <- rbind.fill(lapply(sjc, function(f) {
     data.frame(Filter(Negate(is.null), f), stringsAsFactors=FALSE)
  }))
  rbf
}

#' get data frame for sample
#'
#'
#' @export
getInfoFromRequestSample <- function(sam){
    cn <- function(item){ if(is.null(item)){ NA }else{ item } }
    spikein <- sam$is_spikein
    rq <- sam$request_sample
    rqs <- rq$sample
    pooled <- as.logical(rq$pooled)
    ratio <- rq$ratio
    ownRisk <- as.logical(rqs$own_risk)
    tag <- cn(rqs$tag)
    barcode <- cn(rqs$barcode)
    prep <- rqs$preparation_type
    group <- rqs$group
    exptype <- rqs$exptype
    tag2 <- cn(rqs$secondary_tag)
    multi_id <- cn(rq$multi_id)
    sample_id <- rqs$obj_id
    result <- if(length(sam$check_results) > 0){ cn(sam$check_results[[1]][[1]]$result) }else{ NA }
    data.frame(sampleId=sample_id, multiId=multi_id, isSpikeIn=spikein, barcode=barcode, tag=tag, tag2=tag2, ownRisk=ownRisk, pooled=pooled, ratio=ratio, prep=prep, exptype=exptype, group=group, result=result)
}

#'
#'
#' @export
getSamplesFromFlowcellLane <- function(json){
  sfl <- do.call("rbind", lapply(json$samples, getInfoFromRequestSample  )) 
  sfl
}

#' get samples for flowcell lane
#' only admins can do this!!!  => not completed!!!
#' 
#' @export
getFlowcellLane <- function(flowcell, lane, session){
  s <- NULL
  query <- paste("http://ngs.vbcf.ac.at/forskalle/api/flowcells/", flowcell, "/", lane, sep="")
  tryCatch(
    s <- getURLContent(query, curl=session),
    error=function(e){ cat(paste("error retrieving flowcell/lane info ", flowcell, lane ), file=stderr()) }
  )
  if(is.null(s)){
    return(s)   
  }
  sj <- fromJSON(s) ## its a nested list
  lane_ok <- as.logical(sj$is_ok)
  lane_analyzed <- as.logical(sj$analyzed)
  seqdate <- sj$flowcell$sequencing_date
  seqtype <- sj$flowcell$seqtype
  seqlen <- sj$flowcell$readlen
  seqrap <- sj$flowcell$rapid_mode
  seqpaired <- sj$flowcell$paired
  sams <- getSamplesFromFlowcellLane(sj)
  sams$sequencing_date <- as.Date(seqdate)
  sams$lane_ok <- lane_ok
  sams$lane_analyzed <- lane_analyzed
  sams$seqType <- seqtype
  sams$seqLength <- seqlen
  sams$seqPaired <- seqpaired
  sams$seqRapid <- seqrap
  sams
}


#http://ngs.vbcf.ac.at/forskalle/api/deviceData/request/3901
#http://ngs.vbcf.ac.at/forskalle/api/measurements/request/3901
#' get measurement data for request
#'  
#' @export
getMeasurementsForData <- function(requestId, session){
  s <- NULL
  query <- paste("http://ngs.vbcf.ac.at/forskalle/api/measurements/request/", requestId, sep="")
  tryCatch(
    s <- getURLContent(query, curl=session),
    error=function(e){ cat(paste("error retrieving measurements ", requestId ), file=stderr()) }
  )
  if(is.null(s)){
    return(s)
  }
  sj <- fromJSON(s) ## its a nested list

}

#' extracts from itemList elements with forskalle names and returns a data.frame with 
#' columns named itemNames
#'
#' @export
extractItems <- function(itemList, forskalleNames, itemNames ){
   li <- list()
   #print(str(itemList, 1))
   for( ni in 1:length(forskalleNames)){ 
        nm <- forskalleNames[ni]
        im <- itemNames[ni]
        ex <- itemList[[nm]]
        le <- if(is.null(ex)){ NA }else{ ex }
        #print(paste(nm,im,ex,le, is.null(ex), is.na(ex)))
        li[im] <- le
   }   
   data.frame(li)
}

#' QPCR Extractor
#'
#'
#' @export
extractQPCR <- function(){
   list(form="qPCR",
        type="Data Entry",   
        forskalleNames=c("obj_id", "multi_id", "flag", "efficiency", "form", "user", "size", "date", "control", "kit", "conc", "corrected_conc"),
        itemNames=c("sampleId", "multiId", "flag", "efficiency", "form", "user", "size", "date", "control", "kit", "conc", "corrected_conc")
   )
}

#' cDNA Synthesis Extractor
#' not done
#'
#' @export
extractCDNASynthesis <- function(){
   list(form="cDNA Synthesis",
        type="Data Entry",
        forskalleNames=NA, #c("obj_id", "multi_id", "flag", "efficiency", "form", "user", "size", "date", "control", "kit", "conc", "corrected_conc"),
        itemNames=NA #c("sampleId", "multiId", "flag", "efficiency", "form", "user", "size", "date", "control", "kit", "conc", "corrected_conc")
   )
}

#' extract Size 
#' 
#'
#' @export
extractSize <- function(){
   list(form="Size Analysis",
        type="Data Entry",
        forskalleNames=c("obj_id", "multi_id", "form", "cutout_size", "size", "molarity", "dilution", "conc", "user", "size", "date"),
        itemNames=c("sampleId", "multiId", "form", "cutout size", "size", "molarity", "dilution", "conc", "user", "size", "date")
   )
}


#' extract Preparation
#'
#'
#' @export
extractPreparation <- function(){
   list(form="Preparation",
        type="Data Entry",
        forskalleNames=c("obj_id", "multi_id", "form", "cycles", "kit", "udgase", "date"),
        itemNames=c("sampleId", "multiId", "form", "cycles", "kit", "udgase", "date")
   )
}


#' extract measurement
#'
#'
#' @export
extractMeasurement <- function(measurementList, extractor){
  filtered <- Filter(function(li){ li$form == extractor$form && li$type == extractor$type },  measurementList)
  extracted <- lapply(filtered, extractItems, extractor$forskalleNames, extractor$itemNames)
  df <- do.call("rbind", extracted)
  df[order(df$sampleId, df$date),]
}

#' extract measurements
#'
#'
#' @export
extractMeasurements <- function(measurementList){
  size <- extractMeasurement(measurementList, extractSize())
  preparation <-  extractMeasurement(measurementList, extractPreparation())
  qPCR <- extractMeasurement(measurementList, extractQPCR())
 
  qPCR <- qPCR[!is.na(qPCR$kit),]

  li <- list(size=size, preparation=preparation, qPCR=qPCR)
  li
}



#' get current user details
#' @param session
#'
#' @export
getUserDetails <- function(session){
  s <- NULL
  query <- "http://ngs.vbcf.ac.at/forskalle/api/users/current"
  tryCatch(
    s <- getURLContent(query, curl=session),
    error=function(e){ cat(paste("error retrieving user info"), file=stderr()) }
  )
  if(is.null(s)){
    return(s)
  }
  sj <- fromJSON(s)
  groupName <- sj$primary_group$name
  groups <- data.frame(name=groupName, primary=TRUE)
  list(groups=groups)
}

#' converts NULL to ""
#'
#' @export
NULLtoN <- function(v){
  if(is.null(v)){ "" }else{ v }
}


#' from sample to df
#'
#'  
#' @export
simplifyRunSample <- function(runSample){
   tag1 <- runSample$request_sample$sample$tag
   tag2 <- runSample$request_sample$sample$secondary_tag
   id <- runSample$request_sample$sample$id
   data.frame(id=id,tag1=NULLtoN(tag1),tag2=NULLtoN(tag2))
}

#' from run to df
#'
#' @export
simplifyRunToLaneTags <- function(run){
   flowcell <- run$flowcell_id
   lane <- run$num
   runSamples <- do.call("rbind", lapply(run$samples, simplifyRunSample)) 
   data.frame(flowcell=flowcell, lane=lane, runSamples)
}            


#' get flowcell lanes for logged in user primary group
#' 
#'
#'
#' @export
getRunsForLoggedInGroup <- function(session){
  s <- NULL
  query <- "http://ngs.vbcf.ac.at/forskalle/api/runs/group"
  tryCatch(
    s <- getURLContent(query, curl=session),
    error=function(e){ cat(paste("error retrieving runs info for logged in user "), file=stderr()) }
  )
  if(is.null(s)){
    return(s)
  }
  sj <- fromJSON(s) 
  flowcellLanes <- do.call("rbind", lapply(sj, simplifyRunToLaneTags ))
  flowcellLanes 
}

#' parses bam name into flowcell lane
#'
#'
#' 
#' @export
parseBamName <- function(path){
  p1 <- "([\\w-]*)_(\\d)_(\\d{2}-\\d{2}-\\d{4}).bam"  
  p2 <- "([\\w-]*)_(\\d)_(\\d{4})(\\d{2})(\\d{2})(\\w)_(\\d{8}).bam"
  nm <- basename(path)
  
  exF <- function(nm,p){
    m <- str_match(nm,p)
    if(!is.na(m[1,1])){
       list(flowcell=m[1,2],lane=m[1,3])
    }else{ list() }
  }

  m1 <- exF(nm,p1)
  m2 <- exF(nm,p2)
  if(length(m1) > 0){ m1 }else{ m2 }
}


#' worker method for general users: generate splitfile and command
#' only for first index currently!!
#'
#' on the cluster you have to "module load illumina2bam/1.17"
#'
#' @export
generateSplitFile <- function(bamPath, session){
   fl <- parseBamName(bamPath)
   if(length(fl) > 0){
      outprefix <- paste(fl$flowcell, "_", fl$lane, sep="")
      splitfile <- paste(outprefix, "_barcodes.tab", sep="")
      outdir <- paste(outprefix, "_demultiplexed", sep="")
      metrics <- paste(outdir, "/", outprefix, "_metrics.tab", sep="")      

      fcl <- getRunsForLoggedInGroup(session)      
      tags <- subset(fcl, flowcell==fl$flowcell & lane==fl$lane)
      if(nrow(tags) > 0){
        spl <- data.frame(barcode_sequence=tags$tag1, barcode_name=tags$tag1, library_name=tags$id)
        write.table(spl, splitfile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
        system(paste("mkdir -p ", outdir))    
        cmd <- (paste("BamIndexDecoder OUTPUT_DIR=", outdir, " OUTPUT_PREFIX=", outprefix, " OUTPUT_FORMAT=bam", " INPUT=", bamPath, " BARCODE_FILE=", splitfile, " MAX_RECORDS_IN_RAM=60000000  MAX_MISMATCHES=1 MIN_MISMATCH_DELTA=1 MAX_NO_CALLS=1 CREATE_MD5_FILE=true COMPRESSION_LEVEL=9  METRICS_FILE=", metrics, sep=""))   
        print(cmd)
        cmd
      }else{
        err <- paste("could not find tags for flowcell lane in your primary group runs ", nrow(fcl), fl$flowcell, fl$lane, "\n")
        stop(err)
      }    
   }else{ 
      err <- paste("could not parse flowcell lane from bam ", bamPath, "\n")
      stop(err)
   }  
}






 






