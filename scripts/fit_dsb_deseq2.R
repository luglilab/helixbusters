#!/usr/bin/env Rscript
# Fit one prefiltered feature family; Python owns validation and reporting.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop("Expected counts, metadata, output, design, numerator, denominator")
suppressPackageStartupMessages(library(DESeq2))
counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
metadata <- read.delim(args[2], row.names = 1, check.names = FALSE)
if (!identical(colnames(counts), rownames(metadata))) stop("Sample order differs")
if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) stop("Invalid counts")
storage.mode(counts) <- "integer"
timecourse <- args[4] == "timecourse"
numerators <- if (timecourse) strsplit(args[5], ",", fixed = TRUE)[[1]] else args[5]
metadata$condition <- factor(metadata$group, levels = c(args[6], numerators))
if (anyNA(metadata$condition)) stop("Unexpected conditions")
metadata$donor <- factor(metadata$donor)
formula <- if (args[4] %in% c("paired", "timecourse")) ~ donor + condition else ~ condition
model <- model.matrix(formula, metadata)
if (qr(model)$rank != ncol(model) || nrow(model) <= ncol(model)) stop("Unidentifiable design")
warnings_seen <- character()
withCallingHandlers({
    for (normalization in c("poscounts", "library_total")) {
        dds <- DESeqDataSetFromMatrix(counts, metadata, design = formula)
        if (normalization == "library_total") {
            depth <- metadata$total_molecules
            sizeFactors(dds) <- depth / exp(mean(log(depth)))
        } else {
            dds <- estimateSizeFactors(dds, type = "poscounts")
        }
        # Keep automatic Cook's filtering; do not replace observations with n=3.
        dds <- DESeq(dds, fitType = "parametric", minReplicatesForReplace = Inf,
                     quiet = TRUE, parallel = FALSE)
        if (timecourse) {
            contrasts <- lapply(numerators, function(time) {
                res <- results(dds, contrast = c("condition", time, args[6]),
                               independentFiltering = FALSE, cooksCutoff = TRUE)
                result <- data.frame(feature_id = rownames(res), as.data.frame(res),
                                     dispersion = dispersions(dds),
                                     max_Cooks = apply(assays(dds)[["cooks"]], 1, max),
                                     beta_converged = mcols(dds)$betaConv)
                result$pvalue[is.na(result$beta_converged) | !result$beta_converged] <- NA_real_
                result$padj_within_contrast <- p.adjust(result$pvalue, method = "BH")
                result
            })
            joint <- p.adjust(unlist(lapply(contrasts, function(x) x$pvalue)), method = "BH")
            for (i in seq_along(contrasts)) {
                rows <- seq_len(nrow(contrasts[[i]])) + (i - 1L) * nrow(contrasts[[i]])
                contrasts[[i]]$padj <- joint[rows]
                write.table(contrasts[[i]], file.path(args[3], paste0(normalization, ".", numerators[i], "_vs_", args[6], ".tsv")),
                            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
            }
            # Reuse size factors and dispersions from the all-time-point fit.
            global <- nbinomLRT(dds, reduced = ~ donor, quiet = TRUE)
            res <- results(global, independentFiltering = FALSE, cooksCutoff = TRUE)
            result <- data.frame(feature_id = rownames(res), as.data.frame(res),
                                 dispersion = dispersions(global),
                                 max_Cooks = apply(assays(global)[["cooks"]], 1, max),
                                 full_converged = mcols(global)$fullBetaConv,
                                 reduced_converged = mcols(global)$reducedBetaConv)
            result$pvalue[is.na(result$full_converged) | !result$full_converged |
                          is.na(result$reduced_converged) | !result$reduced_converged] <- NA_real_
            result$padj <- p.adjust(result$pvalue, method = "BH")
            # An omnibus LRT has no single contrast effect size.
            result$log2FoldChange <- NULL
            result$lfcSE <- NULL
            write.table(result, file.path(args[3], paste0(normalization, ".global.tsv")),
                        sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
            write.table(data.frame(sample = colnames(dds), size_factor = sizeFactors(dds)),
                        file.path(args[3], paste0(normalization, ".size_factors.tsv")),
                        sep = "\t", row.names = FALSE, quote = FALSE)
            if (normalization == "poscounts") {
                pdf(file.path(args[3], "dispersion_fit.pdf"))
                plotDispEsts(dds)
                dev.off()
            }
            writeLines(paste(normalization, "dispersion_fit", attr(dispersionFunction(dds), "fitType")),
                       file.path(args[3], paste0(normalization, ".fit.txt")))
            next
        }
        res <- results(dds, contrast = c("condition", args[5], args[6]),
                       independentFiltering = FALSE, cooksCutoff = TRUE)
        result <- data.frame(feature_id = rownames(res), as.data.frame(res),
                             dispersion = dispersions(dds),
                             max_Cooks = apply(assays(dds)[["cooks"]], 1, max),
                             beta_converged = mcols(dds)$betaConv,
                             check.names = FALSE)
        # Failed coefficient convergence is never presented as a valid test.
        result$pvalue[is.na(result$beta_converged) | !result$beta_converged] <- NA_real_
        result$padj <- p.adjust(result$pvalue, method = "BH")
        write.table(result, file.path(args[3], paste0(normalization, ".results.tsv")),
                    sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
        write.table(data.frame(sample = colnames(dds), size_factor = sizeFactors(dds)),
                    file.path(args[3], paste0(normalization, ".size_factors.tsv")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        if (normalization == "poscounts") {
            pdf(file.path(args[3], "dispersion_fit.pdf"))
            plotDispEsts(dds)
            dev.off()
        }
        writeLines(paste(normalization, "dispersion_fit", attr(dispersionFunction(dds), "fitType")),
                   file.path(args[3], paste0(normalization, ".fit.txt")))
    }
}, warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)) })
writeLines(unique(warnings_seen), file.path(args[3], "model_warnings.txt"))
writeLines(c(paste("R", getRversion()), paste("DESeq2", packageVersion("DESeq2")),
             capture.output(sessionInfo())), file.path(args[3], "sessionInfo.txt"))
