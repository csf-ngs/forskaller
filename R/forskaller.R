require("plyr")
require("RCurl")
require("rjson")
require("Hmisc")

## http://ngs.csf.ac.at/forskalle/apidoc

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
#' @param credentials list created with createCredentials
#' @export
startSession <- function(credentials){
  loginurl <- "http://ngs.csf.ac.at/forskalle/api/login"
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
  logouturl <- "http://ngs.csf.ac.at/forskalle/logout" #GET
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
  tryCatch(
   s <- getURLContent(paste("http://ngs.csf.ac.at/forskalle/api/samples/", sampleId, sep=""), curl=session), ## its a string,
   error=function(e){ cat(paste("error retrieving sample info: ", sampleId, "\n", e), file=stderr()) }
  )
  sj <- fromJSON(s) ## its a nested list
  sjl <- removeFromList(sj, "requests")
  logicalColumns <- c("shearing", "add_primer", "own_risk", "fragmented", "stranded")
  logicalColumnsIndex <- names(sjl) %in% logicalColumns 
  sjl[logicalColumnsIndex] <- as.logical(sjl[logicalColumnsIndex])
  # must replace null with NA
  sjln <- nullToNA(sjl)
  
  sampleDF <- as.data.frame(sjln,stringsAsFactors=FALSE)
  sampleDF <- rename(sampleDF, c("preparation_type"="prep"))
  sampleDF$id <- as.integer(sampleDF$id) # is numeric ?? 
  tryCatch(
     m <- getURLContent(paste("http://ngs.csf.ac.at/forskalle/api/measurements/sample/", sampleId, sep=""), curl=session), ## its a string    
     error=function(e){ cat(paste("error retrieving measurement info: ", sampleId, "\n", e), file=stderr()) }
  )
  mj <- fromJSON(m)
  mea <- lapply(mj, measurementToDF)
  if(simplify){
     sm <- lapply(mea, simplifyMeasurement)     
     sa <- simplifySample(sampleDF)
     #sample table is very long, split it somehow in lab/annotation
     #"id"               "preparation_type" "cutout_size"      "shearing"         "fragmented"       "stranded"         "own_risk"         "add_primer"       "exptype"          "organism"         "genotype"         "celltype"         "antibody"         "descr"
     sampleLab <- subset(sa, select=c(id, prep, cutout_size, shearing, fragmented, stranded, own_risk, add_primer))
     sampleAnnot <- subset(sa, select=c(id, exptype, organism, genotype, celltype, antibody, descr))
     list(sample=sa, sampleLab=sampleLab, sampleAnnot=sampleAnnot, measurements=sm)
  }else{
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
   samples <- lapply(sampleIds, getSample, session, TRUE)
   
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
     multi <- getURLContent(paste("http://ngs.csf.ac.at/forskalle/api/multiplexes/", multiId, sep=""), curl=session), ## its a string,
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
 subsetF(preparation, c("obj_id", "batchId", "multi_id", "cycles", "udgase", "cutout_size", "flag", "user", "kit", "method"))  
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
  subs <- subset(sample, select=c(id, tag,  prep, cutout_size, shearing, fragmented, stranded, own_risk, add_primer, exptype, organism, genotype, celltype, antibody, descr))
  within(subs, { genotype=truncateTo(genotype); celltype=truncateTo(celltype); antibody=truncateTo(antibody); descr=truncateTo(descr)})
}


#' simplify RNA Quantification data.frame
#'
#'   resolved conc rin multi_id obj_type                date text notified         change_date  change_user flag       type               form severity    id obj_id         user batchId
#       0  618 8.6     1182   sample 2014-02-20 15:23:32   NA        0 2014-02-21 15:14:57 Carmen Czepe   Ok Data Entry RNA Quantification        0 11278  16864 Carmen Czepe    1376
#'
simplifyRNAQuantification <- function(quantification){
  subsetF(quantification, c("obj_id", "batchId", "multi_id", "conc", "rin", "flag", "user"))
}

#' simplify size analysis data.frame
#' 
#' id severity dilution          form       type        user obj_id notified text                date flag change_user         change_date multi_id obj_type resolved molarity size kit conc method batchId
#' 1 11436        0    -0.51 Size Analysis Data Entry Laura Bayer  16864        0   NA 2014-02-27 14:14:27   Ok Laura Bayer 2014-02-27 15:22:28     1182   sample        0     4.15  270  FA 0.74     HS    1429
simplifySizeAnalysis <- function(sizeanalysis){
  subsetF(sizeanalysis, c("obj_id", "batchId", "multi_id", "dilution", "size", "conc", "flag", "kit", "method", "user")) 
}


#' simplify qpcr data.frame
#'
#'   resolved size efficiency  kit conc multi_id obj_type X2nM_control corrected_conc notified machine text                date flag change_user         change_date severity    id       type form obj_id     user R2 batchId
#'         0  270       96.9 Kapa  2.9     1182   sample         2.24           4.85        0    ours      2014-03-03 11:44:00   Ok    Ru Huang 2014-03-03 11:47:20        0 11478 Data Entry qPCR  16864 Ru Huang  1    1437 
simplifyQPCR <- function(qpcr){
  subsetF(qpcr, c("obj_id", "batchId", "multi_id", "size", "efficiency", "conc", "corrected_conc", "kit", "flag", "user"))
}

#' simplify Quantification (Chip-Seq)
#'
#'
#'      user obj_type volume multi_id notified         change_date bioanalyzer_result obj_id           form resolved change_user  conc       type                date    id total text flag severity bioanalyzer_done batchId
#'     Ru Huang   sample      9     1250        0 2014-03-04 12:52:36                 Ok  17853 Quantification        0    Ru Huang 10010 Data Entry 2014-03-04 11:10:39 11560 90.09   NA   Ok        0                1    1445
simplifyQuantification <- function(quant){
  subsetF(quant, c("obj_id", "batchId", "multi_id", "conc", "flag", "user"))
}



