# Shared phyloseq preparation and random-forest training for the
# microbiota scripts. 01, 02, 03, and export_model.R all use this file
# so the rank tables and models stay the same.

suppressPackageStartupMessages({
  library(phyloseq)
  library(caret)
  library(randomForest)
  library(ggplot2)
})

RANKS <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
RANK_ZH <- c(
  Phylum = "門", Class = "綱", Order = "目",
  Family = "科", Genus = "屬", Species = "種"
)
MAX_FEATURES <- 12L
PREVALENCE_TRY <- c(0.20, 0.10, 0.05)
MODEL_SEED <- 123L
SETTINGS <- list(
  seed = MODEL_SEED,
  max_features = MAX_FEATURES,
  prevalence_try = PREVALENCE_TRY,
  ntree = 500L,
  ranks = RANKS
)

script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)))
  }
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(dirname(normalizePath(ofile, winslash = "/", mustWork = TRUE)))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

is_junk_taxon <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y == "" |
    is.na(x) |
    grepl("unidentified|ambiguous|uncultured|metagenome|unknown|incertae sedis", y, perl = TRUE)
}

safe_colnames <- function(labels) {
  nm <- make.names(labels, unique = TRUE)
  nm[nm == "Group"] <- paste0("taxon_", nm[nm == "Group"])
  nm
}

read_microbiota <- function(root) {
  otu <- as.matrix(read.csv(file.path(root, "otu_table.csv"), row.names = 1, check.names = FALSE))
  storage.mode(otu) <- "double"
  tax <- read.csv(file.path(root, "taxonomy.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  tax[] <- lapply(tax, function(x) trimws(as.character(x)))
  meta <- read.csv(file.path(root, "metadata.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Group" %in% colnames(meta)) stop("metadata.csv needs a Group column")

  common_otu <- intersect(rownames(otu), rownames(tax))
  common_sam <- intersect(colnames(otu), rownames(meta))
  if (length(common_otu) < 2 || length(common_sam) < 10) {
    stop("OTU, taxonomy, and metadata do not share enough samples or taxa")
  }
  otu <- otu[common_otu, common_sam, drop = FALSE]
  tax <- as.matrix(tax[common_otu, , drop = FALSE])
  meta <- meta[common_sam, , drop = FALSE]
  meta$Group <- factor(meta$Group, levels = c("Control", "Parkinson"))
  if (anyNA(meta$Group)) stop("Group contains a label other than Control or Parkinson")

  phy <- phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    tax_table(tax),
    sample_data(meta)
  )
  phy <- prune_taxa(taxa_sums(phy) > 0, phy)
  phy <- prune_samples(sample_sums(phy) > 0, phy)
  phy
}

rank_table <- function(phy, rank) {
  # taxonomyTable's `[` rebuilds a tax_table and refuses a single column.
  # Edit the plain character matrix, then put it back.
  tt <- as(tax_table(phy), "matrix")
  tt[, rank] <- ifelse(is_junk_taxon(tt[, rank]), NA_character_, tt[, rank])
  phy_r <- phyloseq(otu_table(phy), tax_table(tt), sample_data(phy))
  glom <- tax_glom(phy_r, taxrank = rank, NArm = TRUE)
  if (ntaxa(glom) < 2) stop(rank, " has fewer than 2 taxa after tax_glom")

  mat <- as(otu_table(glom), "matrix")
  if (taxa_are_rows(glom)) {
    taxon_ids <- rownames(mat)
    mat <- t(mat)
  } else {
    taxon_ids <- colnames(mat)
  }
  labels <- as.character(as(tax_table(glom), "matrix")[taxon_ids, rank])
  if (anyDuplicated(labels)) stop(rank, " still has duplicated taxon labels after tax_glom")
  colnames(mat) <- labels
  rel <- sweep(mat, 1, rowSums(mat), "/") * 100
  rel[!is.finite(rel)] <- 0
  counts <- mat

  prevalence <- colMeans(counts > 0)
  variance <- apply(rel, 2, stats::var)
  chosen <- NULL
  used_prevalence <- NA_real_
  for (cut in PREVALENCE_TRY) {
    keep <- prevalence >= cut & variance > 0
    if (sum(keep) >= 2) {
      chosen <- names(which(keep))
      used_prevalence <- cut
      break
    }
  }
  if (is.null(chosen)) stop(rank, ": no taxon passed the prevalence filter")
  capped <- FALSE
  if (length(chosen) > MAX_FEATURES) {
    chosen <- names(sort(variance[chosen], decreasing = TRUE))[seq_len(MAX_FEATURES)]
    capped <- TRUE
  } else {
    chosen <- names(sort(variance[chosen], decreasing = TRUE))
  }

  rel <- rel[, chosen, drop = FALSE]
  colnames(rel) <- safe_colnames(chosen)
  df <- as.data.frame(rel, stringsAsFactors = FALSE)
  grp <- as(sample_data(glom), "data.frame")
  df$Group <- factor(grp[rownames(df), "Group"], levels = c("Control", "Parkinson"))

  log_df <- data.frame(
    rank = rank,
    taxon = chosen,
    column = colnames(rel),
    prevalence = as.numeric(prevalence[chosen]),
    variance = as.numeric(variance[chosen]),
    mean_percent = as.numeric(colMeans(rel)),
    kept = TRUE,
    stringsAsFactors = FALSE
  )
  list(
    rank = rank,
    data = df,
    n_glom = ntaxa(glom),
    n_features = ncol(rel),
    prevalence_cutoff = used_prevalence,
    capped = capped,
    feature_log = log_df,
    labels = stats::setNames(chosen, colnames(rel))
  )
}

prepare_microbiota <- function(root) {
  message("Reading OTU, taxonomy, and metadata, then building a phyloseq object...")
  phy <- read_microbiota(root)
  message("phyloseq: ", ntaxa(phy), " OTUs, ", nsamples(phy), " samples")
  tables <- lapply(RANKS, function(rank) {
    message("tax_glom: ", rank)
    rank_table(phy, rank)
  })
  names(tables) <- RANKS
  group <- sample_data(phy)$Group
  list(
    phy_ntaxa = ntaxa(phy),
    phy_nsamples = nsamples(phy),
    group_counts = as.list(table(group)),
    baseline = max(table(group)) / length(group),
    tables = tables
  )
}

train_one_rank <- function(rank_obj) {
  df <- rank_obj$data
  message("Training ", rank_obj$rank, " (", ncol(df) - 1L, " taxa, ", nrow(df), " samples)...")
  set.seed(MODEL_SEED)
  ctrl <- trainControl(method = "cv", number = 10)
  model <- train(Group ~ ., data = df, method = "rf", trControl = ctrl, ntree = 500)
  imp <- varImp(model)
  raw_imp <- importance(model$finalModel)
  if (is.null(dim(raw_imp))) {
    raw <- as.numeric(raw_imp)
    names(raw) <- names(raw_imp)
  } else {
    raw <- as.numeric(raw_imp[, 1])
    names(raw) <- rownames(raw_imp)
  }
  scaled <- as.numeric(imp$importance[, 1])
  names(scaled) <- rownames(imp$importance)
  cols <- colnames(df)[colnames(df) != "Group"]
  raw <- raw[cols]
  scaled <- scaled[cols]
  shifted <- raw - min(raw)
  manual <- shifted / max(shifted) * 100
  if (!isTRUE(all.equal(as.numeric(manual), as.numeric(scaled), tolerance = 1e-8))) {
    stop(rank_obj$rank, ": varImp scaling did not match the caret formula")
  }
  rank_obj$model <- model
  rank_obj$accuracy <- mean(model$resample$Accuracy)
  rank_obj$accuracy_label <- as.character(round(rank_obj$accuracy, 3))
  rank_obj$resample <- model$resample
  rank_obj$best_mtry <- model$bestTune$mtry
  rank_obj$ntree <- model$finalModel$ntree
  rank_obj$importance_raw <- raw
  rank_obj$importance_scaled <- scaled
  rank_obj$tune <- model$results
  message("  CV accuracy ", rank_obj$accuracy_label, "  mtry ", rank_obj$best_mtry)
  rank_obj
}

data_stamp <- function(root) {
  files <- file.path(root, c("otu_table.csv", "taxonomy.csv", "metadata.csv"))
  info <- file.info(files)
  list(files = basename(files), size = as.numeric(info$size), mtime = as.character(info$mtime))
}

build_models <- function(root) {
  prepared <- prepare_microbiota(root)
  trained <- lapply(prepared$tables, train_one_rank)
  acc <- vapply(trained, function(x) x$accuracy, numeric(1))
  best <- names(acc)[which.max(acc)]
  bundle <- list(
    settings = SETTINGS,
    stamp = data_stamp(root),
    phy_ntaxa = prepared$phy_ntaxa,
    phy_nsamples = prepared$phy_nsamples,
    group_counts = prepared$group_counts,
    baseline = prepared$baseline,
    baseline_label = as.character(round(prepared$baseline, 3)),
    accuracy = acc,
    best_rank = best,
    ranks = trained
  )
  saveRDS(bundle, file.path(root, "ml_models.rds"))
  log_df <- do.call(rbind, lapply(trained, function(x) x$feature_log))
  write.csv(log_df, file.path(root, "rank_features.csv"), row.names = FALSE)
  bundle
}

load_models <- function(root, rebuild = FALSE) {
  path <- file.path(root, "ml_models.rds")
  if (!rebuild && file.exists(path)) {
    bundle <- readRDS(path)
    same_settings <- isTRUE(all.equal(bundle$settings, SETTINGS, check.attributes = FALSE))
    same_files <- isTRUE(all.equal(bundle$stamp, data_stamp(root), check.attributes = FALSE))
    if (same_settings && same_files) return(bundle)
    message("Saved models do not match the current data or settings. Retraining.")
  }
  build_models(root)
}

rank_medians <- function(rank_obj, group = NULL) {
  df <- rank_obj$data
  cols <- colnames(df)[colnames(df) != "Group"]
  if (!is.null(group)) df <- df[df$Group == group, , drop = FALSE]
  vapply(df[cols], stats::median, numeric(1))
}

slider_spec <- function(rank_obj) {
  cols <- colnames(rank_obj$data)[colnames(rank_obj$data) != "Group"]
  lapply(cols, function(col) {
    x <- rank_obj$data[[col]]
    mx <- max(x)
    step <- if (mx <= 2) 0.01 else if (mx <= 20) 0.1 else 0.5
    slider_max <- ceiling(mx / step) * step
    if (slider_max <= mx) slider_max <- slider_max + step
    list(
      column = col,
      label = unname(rank_obj$labels[[col]]),
      min = 0,
      max = slider_max,
      step = step,
      value = unname(stats::median(x))
    )
  })
}
