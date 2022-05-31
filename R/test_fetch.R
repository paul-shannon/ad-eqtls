source("fetch.R")
library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_eqtlSummaryCatalog()
   test_eQTLsByLocationAndCategory()
   test_fetch.eqtls.in.chunks()

} # runTests
#----------------------------------------------------------------------------------------------------
test_eqtlSummaryCatalog <- function()
{
   message(sprintf("--- test_eQTLsByLocationAndCategory"))
   brain.studies <- grep("brain", tbl.catalogSummary$unique_id, value=TRUE, ignore.case=TRUE)
   checkTrue(length(brain.studies) > 80)
   checkTrue("GTEx_V8.Brain_Cerebellum" %in% brain.studies)

} # test_eqtlSummaryCatalog
#----------------------------------------------------------------------------------------------------
test_eQTLsByLocationAndCategory <- function()
{
    message(sprintf("--- test_eQTLsByLocationAndCategory"))

    targetGene <- "CLU"
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088

    tbl.cat.geneExpression <- as.data.frame(subset(tbl.catalogSummary, quant_method=="ge"))
    gtex.brain.studies <- grep("gtex_v8.brain", tbl.cat.geneExpression$unique_id,
                               value=TRUE, ignore.case=TRUE)
    cortex.studies <- grep("cortex", gtex.brain.studies, value=TRUE, ignore.case=TRUE)

    chrom <- "8"
    start <- 27603335
    end   <- 27608281
    width <-  1 + end - start

    tbl.1 <- eQTLsByLocationAndStudyID(chrom, start, end, cortex.studies[2], simplify=TRUE)
    tbl.2 <- eQTLsByLocationAndStudyID(chrom, start, end, cortex.studies[2], simplify=TRUE, targetGene="PNOC")

} # test_eQTLsByLocationAndCategory
#----------------------------------------------------------------------------------------------------
test_fetch.eqtls.in.chunks <- function()
{
    message(sprintf("--- test_fetch.eqtls.in.chunks"))

   tbl.small.no.chunks <- fetch.eqtls.in.chunks("chr4", start=23738124, end=23786367,
                                       study="GTEx_V8.Brain_Cerebellar_Hemisphere",
                                       simplify=TRUE, chunk.size=100000)
   checkTrue(nrow(tbl.small.no.chunks) > 200)

   tbl.small.20k.chunks <- fetch.eqtls.in.chunks("chr4", start=23738124, end=23786367,
                                       study="GTEx_V8.Brain_Cerebellar_Hemisphere",
                                       simplify=TRUE, chunk.size=20000)
   checkEquals(nrow(tbl.small.no.chunks), nrow(tbl.small.20k.chunks))

   tbl.small.10k.chunks <- fetch.eqtls.in.chunks("chr4", start=23738124, end=23786367,
                                       study="GTEx_V8.Brain_Cerebellar_Hemisphere",
                                       simplify=TRUE, chunk.size=10000)
   checkEquals(nrow(tbl.small.no.chunks), nrow(tbl.small.10k.chunks))

   if(FALSE){ # big fetch, you probably won't run this routinely
      tbl.big <- fetch.eqtls.in.chunks("chr4", start=22412644, end=25189900,
                              study="GTEx_V8.Brain_Cerebellar_Hemisphere",
                              chunk.size=100000,
                              simplify=TRUE)
      checkTrue(nrow(tbl.big) > 90000)
      } # big fetch

} # test_fetch.eqtls.in.chunks
#----------------------------------------------------------------------------------------------------
