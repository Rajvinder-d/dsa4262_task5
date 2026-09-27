library(bambu)

args <- commandArgs(trailingOnly=TRUE)
se <- readRDS(args[1])
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)

stopifnot(
    all(c("counts", "uniqueCounts", "fullLengthCounts") %in% assayNames(se)),
    "novelTranscript" %in% names(rowData(se)),
    !anyDuplicated(colnames(se))
)
cts <- assay(se, "counts")
uq <- assay(se, "uniqueCounts")
fl <- assay(se, "fullLengthCounts")
novel <- as.logical(rowData(se)$novelTranscript)
exons <- elementNROWS(rowRanges(se))
stopifnot(length(novel) == nrow(se), !anyNA(novel), all(exons >= 1),
          all(is.finite(cts)), all(is.finite(uq)), all(is.finite(fl)),
          all(cts >= 0), all(uq >= 0), all(fl >= 0))
single <- novel & exons == 1
multi <- novel & exons > 1

write_tsv <- function(x, filename) {
    write.table(x, file.path(out, filename), sep="\t", quote=FALSE,
                row.names=FALSE, na="NA")
}

sample_qc <- data.frame(
    sample=colnames(se),
    estimated_counts=round(colSums(cts), 3),
    unique_counts=round(colSums(uq), 3),
    full_length_counts=round(colSums(fl), 3),
    transcripts_ge1=colSums(cts >= 1),
    transcripts_ge10=colSums(cts >= 10),
    transcripts_with_full_length=colSums(fl >= 1),
    novel_single_ge1=colSums(cts[single, , drop=FALSE] >= 1),
    novel_single_ge10=colSums(cts[single, , drop=FALSE] >= 10),
    novel_multi_ge1=colSums(cts[multi, , drop=FALSE] >= 1),
    novel_multi_ge10=colSums(cts[multi, , drop=FALSE] >= 10)
)
write_tsv(sample_qc, "bambu_samples.tsv")

totals <- data.frame(
    metric=c("samples", "total_transcripts", "annotated_transcripts",
             "novel_transcripts", "novel_single_exon_transcripts",
             "novel_multi_exon_transcripts"),
    value=c(ncol(se), nrow(se), sum(!novel), sum(novel), sum(single), sum(multi))
)
write_tsv(totals, "transcript_totals.tsv")

groups <- list(annotated=!novel, novel_single=single, novel_multi=multi)
detected <- rowSums(cts >= 1) > 0
expressed <- rowSums(cts >= 10) > 0
classes <- do.call(rbind, lapply(names(groups), function(group) {
    keep <- groups[[group]]
    data.frame(class=group, total=sum(keep), detected_ge1=sum(keep & detected),
               expressed_ge10=sum(keep & expressed),
               full_length_ge1=sum(keep & rowSums(fl) >= 1),
               unique_ge1=sum(keep & rowSums(uq) >= 1))
}))
write_tsv(classes, "transcripts_by_class.tsv")
writeLines(capture.output(dput(metadata(se)$warnings)),
           file.path(out, "bambu_warnings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(totals, row.names=FALSE)
