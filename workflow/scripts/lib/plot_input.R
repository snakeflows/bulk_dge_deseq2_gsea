# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/lib/plot_input.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-09-21
# Modified: 2026-10-30
# License:  GNU Affero General Public License version 3 or later
# Purpose:  Common MA/volcan plot input code helper R script for the bulk
#           RNA-seq DGE (DESeq2 + GSEA) Snakemake workflow.


###############
# Description #
###############

# Common code used by both `plot_ma.R` and `plot_volcano.R` to optain
# identical input data without code duplication.


###########
# License #
###########

# Common MA/volcan plot input code helper R script for bulk RNA-seq DGE
# (DESeq2 + GSEA) Snakemake workflow.
# Copyright (C) 2026  Marcel Schilling
#
# This file is part of the bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
# workflow.
# The bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake workflow is free software:
# you can redistribute it and/or modify it under the terms of the GNU
# Affero General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

###########
# Imports #
###########

library(readr)
library(dplyr)


#############
# Constants #
#############

GENES_OF_INTEREST_FALLBACK_N <- 10L


##############
# Parameters #
##############

deseq2_results <- snakemake@input[["results"]]
gene_annotation <- snakemake@input[["gene_annotation"]]
contrast_name <- snakemake@params[["contrast_name"]]
lfc_type <- snakemake@params[["lfc_type"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
gene_label <- snakemake@params[["gene_label"]]
genes_of_interest <- snakemake@params[["genes_of_interest"]]
label_non_signif <- snakemake@params[["label_non_signif"]]
color_up   <- snakemake@params[["color_up"]]
color_down <- snakemake@params[["color_down"]]
color_none <- snakemake@params[["color_none"]]

cat("Contrast:", contrast_name, "\n")
cat("LFC type:", lfc_type, "\n")
cat("padj threshold:", padj_threshold, "\n")
cat("Gene label:", gene_label, "\n")
cat("Genes of interest:", genes_of_interest, "\n")
cat("Label non-significant genes:", label_non_signif, "\n")
cat("Colors -- up:", color_up, " down:", color_down, " none:", color_none, "\n\n")

if (is.null(gene_annotation) != is.null(gene_label)) {
  stop("ERROR: Either both or neither of gene annotation and gene label",
       "must be specified!")
}
if (is.null(gene_label)) { gene_label <- "gene_id" }


#########
# Input #
#########

results_df <- read_tsv(deseq2_results, show_col_types = FALSE)

if (gene_label != "gene_id") {
  annotation_df <- read_tsv(gene_annotation, show_col_types = FALSE) %>%
                   select(all_of(c("gene_id", gene_label)))
}


###############
# Gene labels #
###############

if (exists("annotation_df")) {
  if (!(gene_label %in% colnames(annotation_df))) {
    stop("ERROR: Gene annotation does not contain gene label column!")
  }
  results_df <- left_join(results_df, annotation_df, by = "gene_id")
  rm(annotation_df)
  n_missing_labels <- sum(is.na(pull(results_df, gene_label)))
  if (n_missing_labels) {
    cat("WARNING: ", n_missing_labels, " genes have no matching label.\n",
        "Falling back to gene IDs for those genes.\n")
    results_df <- mutate(
      results_df, "{gene_label}" := coalesce(.data[[gene_label]], gene_id)
    )
  }
}


##############
# LFC labels #
##############

lfc_col <- if (lfc_type == "shrunk") "log2FoldChange_shrunk" else "log2FoldChange"
if (!(lfc_col %in% colnames(results_df))) {
  stop("ERROR: column '", lfc_col, "' not found in results table.")
}
if (lfc_type == "shrunk" && all(is.na(results_df[[lfc_col]]))) {
  stop(
    "ERROR: log2FoldChange_shrunk is entirely NA for this ",
    "contrast -- it should not have been requested as a ",
    "shrinkable target. Check the Snakefile's ",
    "SHRINKABLE_CONTRASTS derivation."
  )
}

lfc_label <- if (lfc_type == "shrunk") {
  "Log2 Fold Change (apeglm-shrunk)"
} else {
  "Log2 Fold Change (unshrunk MLE)"
}


#####################
# Genes of interest #
#####################

cat("\nResolving genes of interest...\n")
if (!is.null(genes_of_interest) && length(genes_of_interest) > 0) {
  tested_ids <- results_df$gene_id
  missing_ids <- setdiff(genes_of_interest, tested_ids)
  if (length(missing_ids) > 0) {
    cat(
      "WARNING:", length(missing_ids), "configured genes not found",
      "in this contrast's tested universe; omitting:",
      paste(missing_ids, collapse = ", "), "\n"
    )
  }
  genes_of_interest <- intersect(genes_of_interest, tested_ids)
  if (length(genes_of_interest) == 0) {
    stop(
      "ERROR: none of the configured `genes_of_interest` were found",
      " among the ", nrow(results_df), " gene(s) tested for",
      " contrast '", contrast_name, "'. ",
      "Check the ID(s) in `config[\"genes_of_interest\"]` against",
      " the id column of `results` for this contrast."
    )
  }
  cat(
    "Using", length(genes_of_interest), "configured gene(s) of",
    "interest.\n"
  )
} else {
  cat(
    "No genes of interest configured; falling back to the top",
    GENES_OF_INTEREST_FALLBACK_N, "by padj...\n"
  )
  genes_of_interest <- results_df %>%
    slice_min(padj, n = GENES_OF_INTEREST_FALLBACK_N) %>%
    pull(gene_id)
}


##############
# Annotation #
##############

results_df <- results_df %>%
  mutate(
    lfc_value = .data[[lfc_col]],
    sig = !is.na(padj) & padj < padj_threshold,
    # Note: a significant gene with lfc_value exactly == 0 falls through
    # to "Not Significant" -- effectively impossible with real data, but
    # made explicit rather than left as an implicit case_when fallthrough.
    category = case_when(
      sig & lfc_value > 0 ~ "Up",
      sig & lfc_value < 0 ~ "Down",
      TRUE ~ "Not Significant"
    ),
    category = factor(category, levels = c("Up", "Down", "Not Significant"))
  )

results_df <- results_df %>%
  # Note: empty string rather than `NA` to ensure unlabeled points are
  #  considered for `ggrepel`'s repulsion.
  mutate(label = ifelse((gene_id %in% genes_of_interest) & (
                          label_non_signif | !(category %in% "Not Significant")
                        ), .data[[gene_label]], ""))


##############
# Categories #
##############

category_colors <- c(
  "Up" = color_up,
  "Down" = color_down,
  "Not Significant" = color_none
)

get_category_labels <- function(results_df) {
  category_counts <- table(results_df$category)
  category_labels <- sprintf(
    "%s\n(n=%s)", names(category_counts),
    format(category_counts, big.mark = ",", trim = TRUE)
  )
  names(category_labels) <- names(category_counts)

  category_labels
}


############
# Ordering #
############

order_for_plotting <- function(results_df) {
  # Draw "Not Significant" points first (bottom of the stack) and Up/Down
  # points last (on top), so a dense cloud of gray background points can't
  # visually bury colored hits underneath it. This only reorders the rows
  # handed to ggplot (draw order) -- it has no effect on `category`'s
  # factor levels, which is what controls the legend order.
  results_df %>% arrange(category != "Not Significant")
}
