forskaller
==========

small R API client for our superior in house LIMS

install
-------
library("devtools")
install_github("forskaller", "csf-ngs")


example usage
-------------
```R
library("forskaller")

username <- ""
password <- ""
multiplex <- 0


session <- startSession(createCredentials(username, password))
m <- getMultiplex(multiplex, session)
samples <- getSamples(m$sampleId, session)
samples
samples$measurements$Preparation
samples$measurements$`cDNA Synthesis`
endSession(session)
```





