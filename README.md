forskaller
==========

small R API client for our superior in house LIMS built by Heinz

install
-------
```R
library("devtools")
install_github("forskaller", "csf-ngs")
```

example usage
-------------

##### sample info
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
##### demultiplexing 
[illumina2bam](https://github.com/wtsi-npg/illumina2bam) has to be installed and java -jar /path/to/BamIndexDecoder.jar has to be aliased to BamIndexDecoder

```R
library("forskaller")
username <- ""
password <- ""

session <- startSession(createCredentials(username, password))
cmd <- generateSplitFile(bamPath, session)
print(cmd)
endSession(session)





