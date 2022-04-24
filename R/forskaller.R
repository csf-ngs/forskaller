DEBUG <- FALSE

FSK3ENV <- new.env(parent = emptyenv())
FSK3ENV$API_KEY <- ""
FSK3ENV$DOSTOP <- TRUE


USER_AGENT <- httr::user_agent("http://github.com/hadley/httr")
URLBASE <- "https://ngs.vbcf.ac.at"
FSK3 <- "/forskalle3/api"

#' uses the API-KEY for all interactions
#' either use this or 
#'
#' startSession(createCredentials(username, password))
#' endSession()
#'
#' @param key
#' 
#' @export
useKey <- function(key=NULL){
   if(is.null(key)){
	  key <- yaml::read_yaml("~/.fsk_api.yml")$fsk_api_key
   }
   FSK3ENV$API_KEY <- key
}

#' set dostop on error
#'
#' @export
doStop <- function(st){
  FSK3ENV$DOSTOP <- st
}


#' GET request with optional API-Key and UA
#' 
#' @param path API path 
#' @param query optional query
#'
#' @export
FGET <- function(path, query=NULL, description=NULL){
  apipath <- paste(FSK3, "/", path, sep="")
  apikey <- FSK3ENV$API_KEY 
  hdrs <- if(apikey != ""){
    httr::add_headers(`X-API-Key` = apikey)   
  }else{ NULL }
  if(DEBUG){
    print(paste(apipath, query, "headers:", paste(hdrs, sep=" ", collapse="" )))
  }
  fquery = query[!is.na(query)]
  r <- httr::GET(URLBASE, path=apipath, query=fquery, USER_AGENT, hdrs) 
  stop_if_not_success(r, paste(description, path, query, r$url), FSK3ENV$DOSTOP)
  r
}


##
## https://ngs.vbcf.ac.at/forskalle3/doc/api_intro.md

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
#' the api is a little stateful, but thats behind the back
#' you should call endSession at the end
#' startSession(createCredentials(username, password))
#' 
#' or useKey(API-KEY)
#'
#' @param credentials list created with createCredentials
#'
#' @export
startSession <- function(credentials){
  loginurl <- paste(URLBASE, FSK3, "/login", sep="")
  r <- httr::POST(loginurl, body=credentials, encode="json")
  stop_if_not_success(r, "login")
}


#' ends a forskalle session
#'
#' @export
endSession <- function(){
  logouturl <- "https://ngs.vbcf.ac.at/forskalle3/logout" #GET
  httr::GET(logouturl)
}

#' replace null with NA for as.data.frame
#'
#' @export
nullToNA <- function(li){
  Map(function(v){ if(is.null(v)){ NA }else{ v }}, li)
}

## its a list of individual measurements. The form attribute gives the type ###
## mostly its just one value. changesets is a list => remove, but keep id
## form:
## "Preparation" , "Size Analysis", "qPCR" , "RNA Quantification"
#' converts measurement to df
#'
#' @export
measurementToDF <- function(meas){
  ci <- which(names(meas) == "changeset")
  batchId <- meas$changeset_id
  mf <- meas[-ci]
  mf$batchId <- batchId
  mfn <- nullToNA(mf)
  df <- data.frame(t(rapply(mfn, function(e){ e })), stringsAsFactors=FALSE)
  list(type=mfn$form, data=df)
}



#' removes item by name from list
#'
#' @export
removeFromList <- function(lis, re){
  li <- which(names(lis) == re)
  lis[-li]
}


#' gets sample information from forskalle
#'
#' @param sampleId the sample id 
#' @param session the forskalle session
#' @param sampleTag currently only remove
#' @param should measurement columns be simplified  
#'
#' returns a list of
#' sample = data frame with sample info (1 row)
#' measurments = a list of data frames with measurements 
#'
#' @export
getSample <- function(sampleId, simplify=FALSE, sampleTag="remove"){
  s <- FGET(paste("samples/", sampleId, sep=""))
  sj <- httr::content(s)
  sjl <- sj
  if(sampleTag == "remove"){
     sjl <- removeFromList(sj, "pool_tags")
     sjl <- removeFromList(sjl, "controls")
     sjl <- removeFromList(sjl, "trash_state")
  }
  #logicalColumns <- c("shearing", "add_primer", "own_risk", "fragmented", "stranded")
  logicalColumns <- c("own_risk", "fragmented", "stranded") #was shearing
  logicalColumnsIndex <- names(sjl) %in% logicalColumns 
  nulls <- sapply(sjl[logicalColumnsIndex], is.null)  
  sjl[logicalColumnsIndex][nulls] <- FALSE
  sjl[logicalColumnsIndex] <- as.logical(sjl[logicalColumnsIndex])
  # must replace null with NA
  sjln <- nullToNA(sjl)
  sjln[logicalColumnsIndex][nulls] <- NA
  sampleDF <- as.data.frame(sjln,stringsAsFactors=FALSE)
  sampleDF <- plyr::rename(sampleDF, c("preparation_type"="prep"))
  sampleDF$id <- as.integer(sampleDF$id) # is numeric ?? 
  
  labd <- FGET(paste("samples/", sampleId, "/labdata", sep=""))
  mj <- httr::content(labd)
  mea <- lapply(mj, measurementToDF)
  if(simplify){
     sm <- lapply(mea, simplifyMeasurement)     
     sa <- simplifySample(sampleDF)
     #sample table is very long, split it somehow in lab/annotation
     #"id"               "preparation_type" "cutout_size"      "shearing"         "fragmented"       "stranded"         "own_risk"         "add_primer"       "exptype"          "organism"         "genotype"         "celltype"         "antibody"         "descr"
     sampleLab <- subset(sa, select=c(id, prep, cutout_size, fragmented, stranded, own_risk))#, add_primer))
     sampleAnnot <- subset(sa, select=c(id, exptype, organism, genotype, celltype, antibody, descr))
     list(sample=sa, sampleLab=sampleLab, sampleAnnot=sampleAnnot, measurements=sm)
  } else {
     list(sample=sampleDF, measurements=mea) 
  }
}

#' gets samples as list of data frames
#' 
#' @param sampleIds vector of sample ids
#'  
#' returns a list of data frames
#' samples
#' for each measurement type a data frame
#' all data frames are simplified (as in getSample)
#'
#' @export
getSamples <- function(sampleIds){
   samples <- lapply(sampleIds, getSample, simplify=TRUE)
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
#' @param multiId the id , multiId without M
#'
#' @param session the session
#'
#' @export
getMultiplex <- function(multiId){
   r <- FGET(paste("multiplexes/", multiId, sep=""))
   stop_if_not_success(r, paste("retrieving multiplex", multiId))
   mj <- httr::content(r)
   mjs <- mj$multiplex_samples
   sb <- do.call("rbind", lapply(mjs, function(s){ data.frame(sampleId=s$sample$id, tag=s$sample$adaptor_tag, tag2=s$sample$adaptor_secondary_tag, ratio=s$ratio, stringsAsFactors=FALSE)}))
   sb$multiId <- multiId
   subset(sb,select=c(multiId,sampleId,tag,tag2,ratio)) 
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
    	simplifyUnkown(measurement$data, measurement$type)    
        #stop(paste("unknown measurement ", measurement$type))
    )
    sr <- plyr::rename(simple, c("obj_id"="sampleId", "multi_id"="multiId"))
    list(type=measurement$type, data=sr)
}

#' subsets but is lenient about missing columns
#' missing columns are set to NA
#' @param df the data frame to subset by columns
#' @param cols string vector of columns to subset in desired order
#'
subsetF <- function(df, cols, numcols=NA){
   coli <- match(cols, colnames(df))
 #  if(!is.na(numcols)){
 #     ncoli <- match(numcols, colnames(df))
 #     
 #  }

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
   colnames(alldf) <- gsub("^values.", "", colnames(alldf))
   data.frame(alldf)
}


#' simplifies preparation data.frame
#'
#'   severity    id       type        form         user obj_id cycles udgase text notified cutout_size                date flag  change_user         change_date multi_id obj_type resolved           kit method batchId
#'      0 11344 Data Entry Preparation Carmen Czepe  16864     12      1   NA        0     200-800 2014-02-24 14:00:36   Ok Carmen Czepe 2014-02-26 16:28:39     1182   sample        0 NEB ultra RNA  beads    1391
#'
simplifyPreparation <- function(preparation){
 subsetF(preparation, c("obj_id", "batchId", "multi_id", "values.cycles", "values.method", "values.udgase", "flag", "user", "values.kit", "entry_date"))  
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
    Hmisc::latexTranslate(tru)
  }
  #subs <- subset(sample, select=c(id, tag,  prep, cutout_size, shearing, fragmented, stranded, own_risk, add_primer, exptype, organism, genotype, celltype, antibody, descr))
  subs <- subset(sample, select=c(id, adaptor_tag, adaptor_secondary_tag, prep, cutout_size, fragmented, stranded, own_risk, exptype, organism, genotype, celltype, antibody, scientist, group, description, preparation_kit))
  subs$description <- subs$description ##kweep original
  within(subs, { genotype=truncateTo(genotype); celltype=truncateTo(celltype); antibody=truncateTo(antibody); descr=truncateTo(description)})
}


#' simplify RNA Quantification data.frame
#'
#'   resolved conc rin multi_id obj_type                date text notified         change_date  change_user flag       type               form severity    id obj_id         user batchId
#       0  618 8.6     1182   sample 2014-02-20 15:23:32   NA        0 2014-02-21 15:14:57 Carmen Czepe   Ok Data Entry RNA Quantification        0 11278  16864 Carmen Czepe    1376
#'
simplifyRNAQuantification <- function(quantification){
  subsetF(quantification, c("obj_id", "batchId", "multi_id", "values.conc", "values.kit", "values.rin", "values.total", "values.volume", "flag", "username", "entry_date"))
}

#' simplify size analysis data.frame
#' 
#' id severity dilution          form       type        user obj_id notified text                date flag change_user         change_date multi_id obj_type resolved molarity size kit conc method batchId
#' 1 11436        0    -0.51 Size Analysis Data Entry Laura Bayer  16864        0   NA 2014-02-27 14:14:27   Ok Laura Bayer 2014-02-27 15:22:28     1182   sample        0     4.15  270  FA 0.74     HS    1429
simplifySizeAnalysis <- function(sizeanalysis){
  subsetF(sizeanalysis, c("obj_id", "batchId", "multi_id", "values.conc", "values.dilution", "values.molarity", "values.size", "flag", "values.kit", "method", "username", "entry_date")) 
}

#' simplify unknown, was stop
#'
simplifyUnkown <- function(unknown, mname){
  subsetF(unknown, c("obj_id", "batchId", "multi_id", "username", "entry_date", "unkown measurement", mname))
}

#' simplify qpcr data.frame
#'
#'   resolved size efficiency  kit conc multi_id obj_type X2nM_control corrected_conc notified machine text                date flag change_user         change_date severity    id       type form obj_id     user R2 batchId
#'         0  270       96.9 Kapa  2.9     1182   sample         2.24           4.85        0    ours      2014-03-03 11:44:00   Ok    Ru Huang 2014-03-03 11:47:20        0 11478 Data Entry qPCR  16864 Ru Huang  1    1437 
simplifyQPCR <- function(qpcr){
  subsetF(qpcr, c("obj_id", "batchId", "multi_id", "values.R2", "values.conc", "values.corrected.conc", "values.efficiency", "values.machine", "values.size", "values.kit", "flag", "username", "entry_date"))
}

#' simplify Quantification (Chip-Seq)
#'
#'
#'      user obj_type volume multi_id notified         change_date bioanalyzer_result obj_id           form resolved change_user  conc       type                date    id total text flag severity bioanalyzer_done batchId
#'     Ru Huang   sample      9     1250        0 2014-03-04 12:52:36                 Ok  17853 Quantification        0    Ru Huang 10010 Data Entry 2014-03-04 11:10:39 11560 90.09   NA   Ok        0                1    1445
simplifyQuantification <- function(quant){
  subsetF(quant, c("obj_id", "batchId", "multi_id", "values.bioanalyzer.done", "values.conc", "values.total", "values.volume", "flag", "username", "entry_date"))
}

#' simplify cDNA synthesis (RNA-Seq)
#' 
#'
#'
simplifyCDNASynthesis <- function(cdna){
  subsetF(cdna, c("obj_id", "batchId", "multi_id", "flag", "values.kit", "values.used", "username", "ercc", "entry_date"))
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
createSample <- function(description, comments, group, scientist, celltype="", genotype="", exptype="ChIP-Seq", organism="", tissue_type="", antibody="", status="Aborted", received=today(), fragmented=1, userprep=1, preparation_kit=NULL, preparation_type="none", own_risk=0, barcode="", stranded=0,  secondary_tag=NULL, cutout_size="100-700", tagno=NULL, secondary_tagno=NULL,  ready=NULL, primer="Standard", cutout_size_min=100, cutout_size_max=700, fragment_size=""){
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
runsForSample <- function(sampleId){
  r <- FGET(paste("samples/", sampleId, "/sequencing", sep=""))
  sj <- httr::content(r)
  rbf <- plyr::rbind.fill(lapply(sj, function(f) {
     df <- data.frame(t(rapply(f, function(e){ e })), stringsAsFactors=FALSE)
  }))  
  rbf
}

#'
#'
#' @export
simplifyLaneFromRun <- function(lane){
   fastqccounts <- if(length(lane$lane_checks) > 0 & !is.na(lane$lane_checks[[1]]$fastqcs)){
      lane$lane_checks[[1]]$fastqcs[[1]]$total_count
   } else {
      NA
   }
   samplecount <- lane$sampnum
   exptypes <- sapply(lane$sequenced_samples, function(ss){ ss$request_sample$sample$exptype })   
   exptypesu <- paste(unique(exptypes),collapse=",")
   tibble::tibble(lane=lane$unit_id, exptypes=exptypesu, count=fastqccounts)
}



#' get run by numeric fsk3 id
#'
#' @export
runByFSK3Id <- function(runid){
  r <- FGET(paste("runs/",runid, sep=""))
  sj <- httr::content(r)
  flowcell <- sj$vendor_id
  print(flowcell)
  rows <- do.call("rbind", lapply(sj$run_units, simplifyLaneFromRun))  
  
  rows$flowcell <- sj$vendor_id
  rows$seqdate <- sj$sequencing_date 
  rows
}


#' get samples for own group : available for user
#' if admin, groupName can be from different group
#'  
#'
#' @export
samplesForGroup <- function(groupName, session, since="2017-07-01", admin=FALSE){
  route <- if(!admin){ "samples" }else{ "samples/admin" }
  r <- FGET(route, query=list(filter.received_after=since, filter.group=groupName))
  sj <- httr::content(r)
  rbf <- plyr::rbind.fill(lapply(sj, function(f) {
     df <- data.frame(t(rapply(f, function(e){ e })), stringsAsFactors=FALSE)
  }))
  rbf
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
   plyr::ddply(multiplexDF, .(id), function(s){
      base <- strsplit(s$barcode, "")[[1]]
      position <- seq_along(base)
      laser <- ifelse(base == "A" | base == "C", "red", "green")
      data.frame(id=s$id, base=base, position=position, laser=laser, ratio=s$ratio)
   })
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
#'
#' @export
getFlowcellById <- function(id){
  r <- FGET(paste("runs/illumina/", id, sep=""))
  sj <- httr::content(r)
  rbf <- plyr::rbind.fill(lapply(sj$lanes, function(f) {
     df <- data.frame(t(rapply(f, function(e){ e })), stringsAsFactors=FALSE)
  })) 
  rbf  
}


#' get flowcells  ?from=2014-09-02&to=2015-12-02
#' @param from
#' @param to
#' from <- "2014-09-02"
#'
#' @export
getFlowcells <- function(from="2014-01-02", to=getToday()){
  r <- FGET("runs/illumina", query=list(filter.sequenced_after=from, filter.sequenced_before=to))
  sj <- httr::content(r)  
  sjc <- lapply(sj, removeElementsFromList, c("problems", "lanes", "planned_start", "comments"))    
  rbf <- plyr::rbind.fill(lapply(sjc, function(f) {
     df <- data.frame(t(rapply(f, function(e){ e })), stringsAsFactors=FALSE)
  }))
  rbf
}


#' get requests for ADMIN
#' @param from
#' @param to
#' @param group
#' from <- "2017-01-01"
#'
#' @export
getRequests <- function(group=NA, from="2017-01-01", to=getToday()){
   r <- FGET("requests/admin", query=list(filter.group=group, filter.submitted_after=from, filter.submitted_before=to))
   sj <- httr::content(r)
   requests <- plyr::rbind.fill(lapply(sj, function(f) {
      data.frame(t(rapply(f, function(e){ e })), stringsAsFactors=FALSE)
   }))
   requests
}


#' get sequencings of request 
#' @param id
#'
#' @export
getRequestSequencing <- function(id){
  r <- FGET(paste("requests/",id,"/sequencing",sep=""))
  sj <- httr::content(r)
  sequencings <- plyr::rbind.fill(lapply(sj, function(f) {
     ss1 <- f$sequenced_samples[[1]] #we take the 1. one <<---- is this correct?  
     is_spike <- ss1$is_spikein
     request_lane <- ss1$request_sample$lane
     request_lane_key <- paste(id, "_", request_lane, sep="")
     data.frame(t(rapply(f$run, function(e){ e })),request_id=id, request_lane=request_lane, is_spike=is_spike, sequencing_number=1, request_lane_key, stringsAsFactors=FALSE)
  }))
  if(! is.null(sequencings) && nrow(sequencings) > 1){  
    sequencings <- sequencings[order(sequencings$preparation_date),]
    sequencings$sequencing_number <- 1:nrow(sequencings)
  }
  sequencings
}

#' get estimates of request
#' @param id
#' 
#' @export
getRequestEstimate <- function(id){
   r <- FGET(paste("requests/",id,"/estimate",sep=""))
   sj <- httr::content(r)
   items <- plyr::rbind.fill(lapply(sj$items, function(f) {
     data.frame(request.id=id, request.accepted=NULLtoN(sj$request$accepted), request.submitted=NULLtoN(sj$request$submitted), 
                group=sj$request$group, scientist=sj$request$scientist, 
                category=f$category, code=f$code, count=f$count, description=f$description, price=f$price, total=f$total, 
                cost_assignment=NULLtoN(sj$request$cost_assignment), status=sj$request$status, stringsAsFactors=FALSE)
   }))
   items
}


#' gets request history
#' @param id
#'
#' @export
getRequestHistory <- function(id){
  r <- FGET(paste("requests/",id,"/history",sep=""))
  sj <- httr::content(r)
  hdf <- do.call("rbind", lapply(sj, function(s){ data.frame(t(s), stringsAsFactors=FALSE) }))
  hdf
}

#' date interval to days
#'
#' @export
inDays <- function(t1,t2){
  dt1 <- strftime(t1$date, "%Y-%m-%dT%H:%M:%S%z")
  dt2 <- strftime(t2$date, "%Y-%m-%dT%H:%M:%S%z")
  lubridate::interval(dt1,dt2) / lubridate::ddays(1)
}

#' filter samples
#'
#' filters = list(filter.received_after= , filter.exptype= ) 
#'
#' requires admin 
#'
#'
#' @export
getSamplesByFilter <- function(filters){
  route <- "samples/admin"
  r <- FGET(route, query=filters)
  sj <- httr::content(r)
  rbf <- plyr::rbind.fill(lapply(sj, function(f) {
     df <- data.frame(t(rapply(f, function(e){ e })), stringsAsFactors=FALSE)
  }))
  rbf 
}


###
### I can't go into sample level data easily with the json data
### I would have to parse it, convert to samples, requests_samples +1 for each onlhold -1 for each
### removal of onhold until 0. So now I take 1. onhold until last onhold ends. 
###

#' onhold times in days
#' @export
onholdTimes <- function(request_history_accepted_sequencing){
    onholdStarts <- which(request_history_accepted_sequencing$newstatus == "Onhold")
    onholdEnds <- which(request_history_accepted_sequencing$oldstatus == "Onhold")
    if(length(onholdStarts) > 0 & length(onholdEnds > 0)){ 
       firstStart <- onholdStarts[1]
       lastEnd <- onholdEnds[length(onholdEnds)]
       inDays(request_history_accepted_sequencing[firstStart,]$date, request_history_accepted_sequencing[lastEnd,]$date)       
    } else if(length(onholdStarts) == 0){
       0
    } else if(length(onholdStarts > 0) & length(onholdEnds == 0)) {
       Inf
    }
}

#' accepted to sequence in days
#'
#' @export
fromAcceptedToSequenced <- function(request_history){
  accepted <- which(request_history$newstatus == "Accepted")
  if(length(accepted) > 0){
    sequencing <- which(request_history$newstatus == "Sequencing")   
    if(length(sequencing > 0)){
      request_history_accepted_sequencing <- request_history[accepted[1]:sequencing[1],]
      totaltime <- inDays(request_history_accepted_sequencing$date[1], tail(request_history_accepted_sequencing$date, n=1))
      onhold <- onholdTimes(request_history_accepted_sequencing)
      data.frame(totaltime=totaltime,onhold=onhold,nettime=totaltime-onhold)
    }else{
      data.frame(totaltime=Inf,onhold=NA,nettime=Inf)
    } 
  }else{
    data.frame(totaltime=Inf,onhold=NA,nettime=Inf)
  }
}



prepTypes <- function(request_lane){
   preps <- sapply(request_lane$requests_samples, function(s){ s$sample$preparation_type })
   #pt <- table(preps)
   preps[1]
}

#' get request lanes data
#' @param id
#'
#' @export
getRequest <- function(id){
  r <- FGET(paste("requests/",id,sep=""))
  sj <- httr::content(r)
  reqlanes <- plyr::rbind.fill(lapply(sj$request_lanes, function(f) {
     prep <- prepTypes(f)
     reqL <- extractItems(f, c("status", "share_status", "request_id", "pooled", "multi_id","num"), c("status", "share_status", "request_id", "pooled", "multi_id", "req_lane_num"))
     reqL$request_lane_key <- paste(reqL$request_id, "_", reqL$req_lane_num, sep="")
     reqL$prep <- prep
     reqL
  }))
  reqlanes
}


# TODO: add other columns
# take care not to lose info because factor
getResultCheckData <- function(resultCheck){
   rcd <- data.frame(resultCheck[[1]], stringsAsFactors=FALSE)
   srcd <- rcd[c("total", "countq30", "dupl")]
   as.data.frame(lapply(srcd, as.numeric))
   #srcd <- subsetF(rcd, c("total", "countq30", "dupl", "bloomcat"))
   #as.data.frame(lapply( ,as.numeric))
}

#' get data frame for sample
#'
#'
#' @export
getInfoFromRequestSample <- function(sam, withResults=FALSE){
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
    sampleInfo <- data.frame(sampleId=sample_id, multiId=multi_id, isSpikeIn=spikein, barcode=barcode, tag=tag, tag2=tag2, ownRisk=ownRisk, pooled=pooled, ratio=ratio, prep=prep, exptype=exptype, group=group, result=result)

    if(withResults){
      rc <- getResultCheckData(sam$check_results[[1]])
      sr <- cbind(sampleInfo,rc)
      sr
    }else{
      sampleInfo
    }

}

#'
#'
#' @export
getSamplesFromFlowcellLane <- function(json, withResult){
  sampleInfo <-  lapply(json$samples, getInfoFromRequestSample, withResult)
  sfl <- do.call("rbind", sampleInfo)
  sfl
}



perLane <- function(lane){
    num <- as.integer(lane$num)
    lane_ok <- as.logical(lane$is_ok)
    lane_analyzed <- as.logical(lane$analyzed)
    primer <- lane$primer
    withData <- FALSE
    total <- NA
    countq30 <- NA
    if(length(lane$unsplit_checks) > 0){
        total <- as.integer(lane$unsplit_checks[[1]]$total)
        countq30 <- as.integer(lane$unsplit_checks[[1]]$countq30)
        withData <- TRUE
    }
    data.frame(lane=num, lane_ok=lane_ok, lane_analyzed=lane_analyzed, primer=primer, total=total, q30=countq30, withData=withData)            
}


getLaneChecks <- function(lane){
  nr <- lane$unit_id 
  if(length(lane$lane_checks) > 0){
    fqcs <- lane$lane_checks[[1]]$fastqcs
    df <- data.frame(t(rapply(f, function(e){ e })), stringsAsFactors=FALSE)
    df$nr <- nr
  }else{
    df <- data.frame(nr=nr)
  } 
  df
}

#' get stats for flowcell lane without samples
#' only admins can do this!!!
#' 
#' @export
getFlowcellStats <- function(flowcell){
  r <- FGET(paste("runs/illumina/", flowcell, sep="")) 
  fc <- httr::content(r)
  lanes <- fc$lanes
  rbf <- plyr::rbind.fill(lapply(lanes, function(f) {
     getLaneChecks(f)
  }))
  rbf$seqdate <- fc$sequencing_date
  rbf
}

# for each sample:
getLaneSample <- function(sequencedSample){
   id <- sequencedSample$request_sample$sample$id
   multi_id <- sequencedSample$multi_id
   rs <- sequencedSample$request_sample
   rl <- rs$request_lane
   group <- rl$request$group
   groupId <- rl$request$group_id
   demultiplexing <- rl$request$demultiplexing
   share_status <- rl$share_status
   share_required_ratio <- rl$share_required_ratio
   request_id <- rs$request_id   
   is_spikein <- sequencedSample$is_spikein
   tibble::tibble(id=id,demultiplexing=demultiplexing, is_spikein=is_spikein, multi_id=NULLtoNA(multi_id),request_id=request_id, share_status=NULLtoNA(share_status), share_required_ratio=NULLtoNA(share_required_ratio), group=group, group_id=groupId) 
}



getLaneSamples <- function(lane){
  ss <- lane$sequenced_samples
  samples <-  plyr::ldply(ss,getLaneSample) 
  samples$lane <- lane$num
  samples$lane_status <- lane$status
  samples
}

#' get lanes for flowcell
#' 
#' only admins
#' @export
getFlowcellLanes <- function(flowcell){
  print(flowcell)
  r <- FGET(paste("runs/illumina/", flowcell, sep=""))
  fc <- httr::content(r)
  lanes <- fc$lanes
  rbf <- plyr::rbind.fill(lapply(lanes, function(f) {
     getLaneSamples(f)
  }))
  rbf$seqdate <- fc$sequencing_date
  rbf$flowcell <- flowcell
  rbf
}


#http://ngs.vbcf.ac.at/forskalle/api/deviceData/request/3901
#http://ngs.vbcf.ac.at/forskalle/api/measurements/request/3901
#' get measurement data for request
#'  
#' @export
getMeasurementsForRequest <- function(requestId, session){
  stop("not done")
  
  #r <- FGET(paste("measurements/request", requestId, sep=""))
  #fc <- httr::content(r)
  #rbf <- plyr::rbind.fill(lapply(lanes, function(f) {
     
  #}))
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

#' converts NULL to NA
#'
#' @export
NULLtoNA <- function(v){
  if(is.null(v)){ NA }else{ v }
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

#' converts a list to a data frame, setting NULL to NA, only flat lists
#'
#' @export
listToDF <- function(li){
  liNA <- Map(function(le){ if(is.null(le)){ le <- NA }else{ le }}, li)
  as.data.frame(liNA, stringsAsFactors=FALSE)
}

#' random => ""
#'
#' @export
removeRandom <- function(column){
  column[column == "random"] <- ""
  column
}


#' get barcodes for flowcell lane
#' 
#' @param flowcell the flowcell
#' @param lane the lane
#'
#' @export 
getBarcodes <- function(flowcell, lane){
   r <- FGET(paste("runs/illumina/", flowcell, "/", lane, "/barcodes", sep=""))
   mj <- httr::content(r)
   fb <- do.call("rbind", lapply(mj, listToDF))
   fb$adaptor_tag <- removeRandom(fb$adaptor_tag)
   fb$adaptor_secondary_tag <- removeRandom(fb$adaptor_secondary_tag)  
   fb
}

#' get barcodes file for flowcell lane 
#'
#' @param flowcell the flowcell
#' @param lane the lane
#' @param outpath output path for tab delimited file
#'
#' @export 
writeBarcodesFile <- function(flowcell, lane, outpath){
   tab <- getBarcodes(flowcell, lane)
   writeBarcodesToFile(tab, outpath)
}
 
#' writes barcodes to dual index file
#'
#' @param outpath
#' @param barcodes
#' @export
writeBarcodesToFile <- function(barcodes, outpath){
   cols <- c("adaptor_tag", "adaptor_secondary_tag", "sample_id")
   tabs <- subset(barcodes, select=cols)
   na2 <- is.na(tabs$adaptor_secondary_tag)
   tabs$adaptor_secondary_tag[na2] <- ""
   tabs$outname <- paste(tabs$sample_id, "_", tabs$adaptor_tag, tabs$adaptor_secondary_tag, sep="")
   write.table(tabs, outpath, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}


#' stops if response is not successs
#'
#' @param response (the http respons)
#' @param message the optional error message
#' 
#' @export
stop_if_not_success <- function(response, message="", dostop=TRUE){
   if(httr::http_status(response)$category != "Success"){
     error <- as.character(httr::http_status(response))
     if(message != ""){
        error <- paste(response[1], error, message, sep=" ", collapse="\n")
     }
     write(error, stderr())
     if(DEBUG){
       write(paste("headers:",as.data.frame(response$request$headers)), stderr())
     }
     if(dostop){
        stop(error)
     }
   }
}


 
#' get plates 
#' 
#' @return //list of plates
#'
#' @export
get_plates_for_request <- function(request_id){
	r <- FGET(paste("requests/", request_id, "/plates", sep=""))
    httr::content(r)
}	

#' plate cell to tibble
#'
#' @export
plate_cell_to_tibble <- function(cell){
	sid <- if(is.null(cell$sample)){ NA }else{ as.character(cell$sample$id) }
	tibble::tibble(row=cell$row, column=cell$column, sample_id=sid)
}

#' plate row to tibble
#'
#' @export
plate_row_to_tibble <- function(row){
	cells <- lapply(row$cells, plate_cell_to_tibble) 
 	do.call("rbind", cells)
} 


#' plate to tibble
#' sample row column plate
#'
#' @export
plate_layout_to_tibble <- function(plate){
	rows <- lapply(plate$rows, plate_row_to_tibble)   	
	tr <- do.call("rbind", rows)
	tr$plate_type <- plate$type
	tr
}



