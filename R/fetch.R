library(catalogueR)   # make sure you build the dev branch:
library(EnsDb.Hsapiens.v79)
library(plyr)  # for rbind.fill
tbl.catalogSummary <- get(data(meta, package="catalogueR"))
#----------------------------------------------------------------------------------------------------
eQTLsByLocationAndStudyID <- function(chrom, start, end, studyIDs, targetGene=NA, simplify=FALSE)
{
    tbls <- list()
    for(id in studyIDs){
        message(sprintf("--- fetching %s (ge)", id))
        tryCatch({
            #browser()
            suppressWarnings({tbl <- eQTL_Catalogue.fetch(unique_id=id,
                                                          quant_method="ge",
                                                          method="REST",
                                                          chrom = sub("chr", "", chrom),
                                                          bp_lower=start,
                                                          bp_upper=end,
                                                          verbose=TRUE)})
            printf("%s %d-%d, %d", id, start, end, nrow(tbl))
            if(nrow(tbl) > 0){
                tbl$id <- id
                tbls[[id]] <- tbl
            }
        },
        error = function(e){
            message(sprintf("eQTL_Catalogue.fetch failed on study %s", id))
            print(e)
        })
    } # for id
    tbl.out <- do.call(rbind.fill, tbls)
    rownames(tbl.out) <- NULL
    if(is.null(tbl.out))
        return(data.frame())
    new.order <- order(tbl.out$pvalue.QTL, decreasing=FALSE)
    tbl.out <- as.data.frame(tbl.out[new.order,])
    coi <- c("rsid.QTL", "pvalue.QTL", "gene_id.QTL", "an.QTL", "beta.QTL", "id")
    if(simplify){
        tbl.out <- tbl.out[, coi]
        colnames(tbl.out) <- c("rsid", "pvalue", "gene", "total.alleles", "beta", "id")
        }
    if(!is.na(targetGene)){
        ensg <- as.character(mapIds(EnsDb.Hsapiens.v79, targetGene, "GENEID", "SYMBOL"))
        if(!ensg %in% tbl.out$gene){   # these eQTLs are not for the target gene
            message(sprintf("ADv$geteQTLsByLocationAndStudyID, %d eqtls found, but none for %s",
                            nrow(tbl.out), targetGene))
            return(data.frame())
        }
        tbl.out <- subset(tbl.out, gene==ensg)
        tbl.out$gene <- targetGene
        message(sprintf("--- %d variants for %s, corrected from %s", nrow(tbl.out), targetGene, ensg))
    } else {
        map <- mapIds(EnsDb.Hsapiens.v79, tbl.out$gene, "SYMBOL", "GENEID")
        tbl.map <- data.frame(ensg=names(map), symbol=as.character(map), stringsAsFactors=FALSE)
        na.indices <- which(is.na(tbl.map$symbol))
        length(na.indices)
        tbl.map$symbol[na.indices] <- tbl.map$ensg[na.indices]
        tbl.out$gene <- tbl.map$symbol
    }
    rownames(tbl.out) <- NULL
    invisible(as.data.frame(tbl.out))

} # eQTLsByLocationAndCategory
#----------------------------------------------------------------------------------------------------
fetch.eqtls.in.chunks <- function(chrom, start, end, study, simplify, chunk.size)
{
  roi.width <- 1 + end - start

  if(roi.width <= chunk.size){
     message(sprintf("--- pgc1a fetch.eqtls, just one chunk"))
     tbl <- eQTLsByLocationAndStudyID(chrom, start, end, study, simplify=simplify)

  } else {
     boundaries.needed <- 1 + (roi.width %/% chunk.size)
     starts <- as.integer(seq(from=start, to=end, length.out=boundaries.needed))
     ends <- starts[-1]
     starts <- starts[-(length(starts))]
     tbls <- list()
     intervals <- length(starts)
     message(sprintf("==== pgc1a fetch.eqtls, %d chunks", intervals))
     for(i in seq_len(intervals)){
        message(sprintf("--- fetching chunk %2d/%d for %s", i, intervals, study))
        tbl.chunk <- eQTLsByLocationAndStudyID(chrom,
                                               as.integer(starts[i]),
                                               as.integer(ends[i]),
                                               study,
                                               targetGene=NA,
                                               simplify=simplify)
        tbls[[i]] <- tbl.chunk
        } # for i
     tbl <- do.call(rbind, tbls)
     } # else

  invisible(tbl)

} # fetch.eqtls.in.chunks
#----------------------------------------------------------------------------------------------------
