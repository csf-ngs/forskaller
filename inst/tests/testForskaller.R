context("forskaller tests")

T <- TRUE
F <- FALSE


test_that("flowcell parsing", {
    path1 <- "/some/junk/path/C9B6NANXX_8_20-03-2016.bam"
    path2 <- "/other/junk/path/C9B7LANXX_1_20160322B_20160323.bam"
    path3 <- "C9B6NANXX_3_20-03-2016.bam"
    path4 <- "/one/myseq/name/000000000-AHJUJ_1_01-10-2015.bam"
    wrong <- "/one/myseq/name/000000000-AHJUJ_X_01-10-2015.bam"
    expect_equal(parseBamName(path1), list(flowcell="C9B6NANXX", lane="8"))
    expect_equal(parseBamName(path2), list(flowcell="C9B7LANXX", lane="1"))
    expect_equal(parseBamName(path3), list(flowcell="C9B6NANXX", lane="3"))
    expect_equal(parseBamName(path4), list(flowcell="000000000-AHJUJ", lane="1"))
    expect_equal(parseBamName(wrong), list())
})
