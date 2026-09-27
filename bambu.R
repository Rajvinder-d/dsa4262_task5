library(bambu)

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) >= 5L, args[3] %in% c("annotated", "unannotated"))
fa.file <- args[1]
annotation.mode <- args[3]
annotations <- if (annotation.mode == "annotated") prepareAnnotations(args[2]) else NULL
cores <- as.integer(args[4])
samples.bam <- args[-(1:4)]

settings <- c(paste0("annotation_mode=", annotation.mode),
              paste0("annotations=", if (is.null(annotations)) "NULL" else args[2]),
              "NDR=1", "min.readFractionByEqClass=0.2",
              "min.txScore.singleExon=0", paste0("ncore=", cores))
writeLines(settings, "bambu/run_settings.txt")
cat(paste(settings, collapse="\n"), "\n")

writeLines(capture.output(sessionInfo()), "bambu/sessionInfo.txt")
writeLines(basename(samples.bam), "bambu/samples.txt")

set.seed(1)
se <- bambu(
    reads=samples.bam,
    annotations=annotations,
    genome=fa.file,
    opt.discovery=list(
        min.txScore.singleExon=0
    ),
    ncore=cores,
    lowMemory=FALSE,
    verbose=TRUE
)

saveRDS(se, "bambu/se.rds")
writeBambuOutput(se, path="bambu")
# Stable, explicitly named transcript read-count matrix for the assignment.
write.table(data.frame(transcript_id=rownames(se), assay(se, "counts"),
                       check.names=FALSE),
            "bambu/read_counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)
show(se)
