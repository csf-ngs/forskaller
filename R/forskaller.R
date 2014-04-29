require("plyr")
require("RCurl")
require("rjson")
require("Hmisc")

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
  loginResult <- postForm(loginurl, .params = credentials, curl=curl, .checkParams=FALSE )
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
  s <- getURLContent(paste("http://ngs.csf.ac.at/forskalle/api/samples/", sampleId, sep=""), curl=session) ## its a string
  sj <- fromJSON(s) ## its a nested list
  sjl <- removeFromList(sj, "requests")
  logicalColumns <- c("shearing", "add_primer", "own_risk", "fragmented", "stranded")
  logicalColumnsIndex <- names(sjl) %in% logicalColumns 
  sjl[logicalColumnsIndex] <- as.logical(sjl[logicalColumnsIndex])
  # must replace null with NA
  sjln <- nullToNA(sjl)
  
  sampleDF <- as.data.frame(sjln,stringsAsFactors=FALSE)
  
  m <- getURLContent(paste("http://ngs.csf.ac.at/forskalle/api/measurements/sample/", sampleId, sep=""), curl=session) ## its a string    
  mj <- fromJSON(m)
  mea <- lapply(mj, measurementToDF)
  if(simplify){
     sm <- lapply(mea, simplifyMeasurement)     
     sa <- simplifySample(sampleDF)
     list(sample=sa, measurements=sm)
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
   samplesDF <- do.call("rbind", lapply(samples, function(s){ s$sa }))
   measurementTypesM <- sapply(samples, function(s){  
      sapply(s$measurements, function(m){ m$type })
   })
   measurementTypes <- unique(as.vector(measurementTypesM))

   measurements <- lapply(measurementTypes, function(type){
        meas <- Map(function(s){ s$measurements }, samples)
        ms <- Filter(function(m){ m$type == type }, unlist(meas, recursive=FALSE))
        dat <- do.call("rbind", Map(function(m){ m$data }, ms))
        dato <- dat[order(dat$sampleId),] 
        dato
   })
   names(measurements) <- measurementTypes    
   sampleso <- samplesDF[order(samplesDF$id),]
   list(samples=sampleso, measurements=measurements)
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
        'RNA Quantification'=simplifyRNAQuantification(measurement$data)
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
 subsetF(preparation, c("id", "obj_id", "batchId", "multi_id", "cycles", "udgase", "cutout_size", "flag", "user", "kit", "method"))  
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
  subs <- subset(sample, select=c(id, tag, preparation_type, cutout_size, shearing, fragmented, stranded, own_risk, add_primer, exptype, organism, genotype, celltype, antibody, descr))
  within(subs, { genotype=truncateTo(genotype); celltype=truncateTo(celltype); antibody=truncateTo(antibody); descr=truncateTo(descr)})
}


#' simplify RNA Quantification data.frame
#'
#'   resolved conc rin multi_id obj_type                date text notified         change_date  change_user flag       type               form severity    id obj_id         user batchId
#       0  618 8.6     1182   sample 2014-02-20 15:23:32   NA        0 2014-02-21 15:14:57 Carmen Czepe   Ok Data Entry RNA Quantification        0 11278  16864 Carmen Czepe    1376
#'
simplifyRNAQuantification <- function(quantification){
  subsetF(quantification, c("id", "obj_id", "batchId", "multi_id", "conc", "rin", "flag", "user"))
}

#' simplify size analysis data.frame
#' 
#' id severity dilution          form       type        user obj_id notified text                date flag change_user         change_date multi_id obj_type resolved molarity size kit conc method batchId
#' 1 11436        0    -0.51 Size Analysis Data Entry Laura Bayer  16864        0   NA 2014-02-27 14:14:27   Ok Laura Bayer 2014-02-27 15:22:28     1182   sample        0     4.15  270  FA 0.74     HS    1429
simplifySizeAnalysis <- function(sizeanalysis){
  subsetF(sizeanalysis, c("id", "obj_id", "batchId", "multi_id", "dilution", "size", "conc", "flag", "kit", "method", "user")) 
}


#' simplify qpcr data.frame
#'
#'   resolved size efficiency  kit conc multi_id obj_type X2nM_control corrected_conc notified machine text                date flag change_user         change_date severity    id       type form obj_id     user R2 batchId
#'         0  270       96.9 Kapa  2.9     1182   sample         2.24           4.85        0    ours      2014-03-03 11:44:00   Ok    Ru Huang 2014-03-03 11:47:20        0 11478 Data Entry qPCR  16864 Ru Huang  1    1437 
simplifyQPCR <- function(qpcr){
  subsetF(qpcr, c("id", "obj_id", "batchId", "multi_id", "size", "efficiency", "conc", "corrected_conc", "kit", "flag", "user"))
}





