# Export every rank's random forest into index.html.
# The browser walks the same trees as randomForest::predict.

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
root <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(root, "microbiota_common.R"))
suppressPackageStartupMessages({
  library(jsonlite)
  library(V8)
})

f64_b64 <- function(x) {
  if (!length(x)) return("")
  con <- rawConnection(raw(0), "w")
  writeBin(as.double(x), con, size = 8, endian = "little")
  val <- rawConnectionValue(con)
  close(con)
  gsub("\\s+", "", jsonlite::base64_enc(val))
}

extract_trees <- function(rf, nclass) {
  lapply(seq_len(rf$ntree), function(k) {
    tr <- getTree(rf, k = k, labelVar = FALSE)
    status_terminal <- tr[, 5] == -1
    if (!identical(as.logical(status_terminal), as.logical(tr[, 3] == 0))) {
      stop("Tree ", k, ": terminal flag mismatch")
    }
    pred <- as.integer(tr[, 6])
    terminal <- which(status_terminal)
    if (!length(terminal) || any(pred[terminal] < 1L | pred[terminal] > nclass)) {
      stop("Tree ", k, ": leaf class out of range")
    }
    list(
      L1 = as.integer(tr[, 1]),
      R1 = as.integer(tr[, 2]),
      L = as.integer(ifelse(tr[, 1] == 0, -1L, tr[, 1] - 1L)),
      R = as.integer(ifelse(tr[, 2] == 0, -1L, tr[, 2] - 1L)),
      V = as.integer(tr[, 3]),
      S = as.numeric(tr[, 4]),
      C = pred
    )
  })
}

walk_one <- function(tr, x) {
  node <- 1L
  guard <- 0L
  while (tr$V[node] != 0L) {
    guard <- guard + 1L
    if (guard > length(tr$V)) stop("walker did not reach a leaf")
    if (x[[tr$V[node]]] <= tr$S[node]) node <- tr$L1[node] else node <- tr$R1[node]
  }
  tr$C[node]
}

bundle <- load_models(root)
message("Exporting ", length(RANKS), " forests...")

export_rank <- function(rank) {
  obj <- bundle$ranks[[rank]]
  rf <- obj$model$finalModel
  classes <- as.character(rf$classes)
  vnames <- rownames(rf$importance)
  if (is.null(vnames)) vnames <- colnames(obj$data)[colnames(obj$data) != "Group"]
  trees <- extract_trees(rf, length(classes))
  check_df <- head(obj$data[, vnames, drop = FALSE], 5)
  pa <- predict(rf, check_df, predict.all = TRUE)
  indiv <- matrix(as.character(pa$individual), nrow = nrow(check_df), ncol = rf$ntree)
  for (i in seq_len(nrow(check_df))) {
    x <- as.numeric(check_df[i, vnames])
    for (k in seq_len(rf$ntree)) {
      got <- classes[walk_one(trees[[k]], x)]
      if (!identical(got, indiv[i, k])) {
        stop(rank, " tree walk mismatch at row ", i, " tree ", k)
      }
    }
  }
  spec <- slider_spec(obj)
  spec <- spec[match(vnames, vapply(spec, function(f) f$column, character(1)))]
  features <- lapply(spec, function(f) {
    list(
      column = unbox(f$column),
      label = unbox(f$label),
      min = unbox(f$min),
      max = unbox(f$max),
      step = unbox(f$step),
      value = unbox(f$value)
    )
  })
  raw <- unname(obj$importance_raw[vnames])
  scaled <- unname(obj$importance_scaled[vnames])
  rs <- obj$resample[order(obj$resample$Resample), , drop = FALSE]
  folds <- lapply(seq_len(nrow(rs)), function(i) {
    list(name = unbox(as.character(rs$Resample[i])), Accuracy = unbox(rs$Accuracy[i]))
  })
  tune <- lapply(seq_len(nrow(obj$tune)), function(i) {
    list(
      mtry = unbox(as.integer(obj$tune$mtry[i])),
      Accuracy = unbox(obj$tune$Accuracy[i]),
      selected = unbox(as.integer(obj$tune$mtry[i]) == as.integer(obj$best_mtry))
    )
  })
  med_all <- rank_medians(obj)
  med_c <- rank_medians(obj, "Control")
  med_p <- rank_medians(obj, "Parkinson")
  presets <- list(
    list(name = unbox("全體中位數"), x = I(unname(med_all[vnames]))),
    list(name = unbox("Control 中位數"), x = I(unname(med_c[vnames]))),
    list(name = unbox("Parkinson 中位數"), x = I(unname(med_p[vnames])))
  )
  train_pred <- predict(obj$model, obj$data)
  tab <- table(
    Predicted = factor(train_pred, levels = classes),
    Actual = factor(obj$data$Group, levels = classes)
  )
  list(
    meta = list(
      rank = unbox(rank),
      zh = unbox(unname(RANK_ZH[[rank]])),
      classes = I(unname(classes)),
      cutoff = I(unname(as.numeric(rf$forest$cutoff))),
      features = features,
      ntree = unbox(as.integer(rf$ntree)),
      mtry = unbox(as.integer(rf$mtry)),
      accuracy = unbox(obj$accuracy),
      accuracyLabel = unbox(obj$accuracy_label),
      accuracyText = unbox(format(obj$accuracy, digits = 16)),
      nGlom = unbox(obj$n_glom),
      nFeatures = unbox(obj$n_features),
      prevalence = unbox(obj$prevalence_cutoff),
      capped = unbox(obj$capped),
      folds = folds,
      tune = tune,
      importance = list(
        features = I(unname(vnames)),
        labels = I(unname(obj$labels[vnames])),
        raw = unbox(f64_b64(raw)),
        scaled = I(unname(as.numeric(scaled)))
      ),
      presets = presets,
      confusion = list(
        classes = I(unname(classes)),
        matrix = lapply(seq_along(classes), function(i) I(unname(as.integer(tab[i, ]))))
      ),
      trainAccuracyLabel = unbox(as.character(round(mean(train_pred == obj$data$Group), 3))),
      featureLog = lapply(seq_len(nrow(obj$feature_log)), function(i) {
        list(
          taxon = unbox(obj$feature_log$taxon[i]),
          prevalence = unbox(obj$feature_log$prevalence[i]),
          mean = unbox(obj$feature_log$mean_percent[i])
        )
      })
    ),
    trees = lapply(trees, function(tr) {
      list(
        L = I(unname(tr$L)),
        R = I(unname(tr$R)),
        V = I(unname(tr$V)),
        S = unbox(f64_b64(tr$S)),
        C = I(unname(tr$C))
      )
    }),
    vnames = vnames,
    classes = classes,
    data = obj$data,
    model = obj$model,
    rf = rf
  )
}

pieces <- lapply(RANKS, export_rank)
names(pieces) <- RANKS

model_obj <- list(
  rankOrder = I(RANKS),
  rankZh = lapply(RANK_ZH[RANKS], unbox),
  baseline = unbox(bundle$baseline),
  baselineLabel = unbox(bundle$baseline_label),
  nSamples = unbox(bundle$phy_nsamples),
  nOtu = unbox(bundle$phy_ntaxa),
  nControl = unbox(unname(bundle$group_counts$Control)),
  nParkinson = unbox(unname(bundle$group_counts$Parkinson)),
  bestRank = unbox(bundle$best_rank),
  maxFeatures = unbox(MAX_FEATURES),
  versions = list(
    r = unbox(R.version.string),
    caret = unbox(as.character(packageVersion("caret"))),
    randomForest = unbox(as.character(packageVersion("randomForest"))),
    phyloseq = unbox(as.character(packageVersion("phyloseq"))),
    rng = unbox(paste(RNGkind(), collapse = " / ")),
    seed = unbox(MODEL_SEED)
  ),
  generated = unbox(format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  ranks = lapply(pieces, function(p) {
    out <- p$meta
    out$trees <- p$trees
    out
  })
)

message("Checking JavaScript against R...")
predict_src <- paste(readLines(file.path(root, "predict.js"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
ctx <- v8()
invisible(ctx$eval(predict_src))
verify_ranks <- list()

for (rank in RANKS) {
  p <- pieces[[rank]]
  set.seed(456)
  spec <- slider_spec(bundle$ranks[[rank]])
  spec <- spec[match(p$vnames, vapply(spec, function(f) f$column, character(1)))]
  rand <- as.data.frame(lapply(spec, function(f) runif(24, f$min, f$max)))
  colnames(rand) <- p$vnames
  preset_mat <- do.call(rbind, lapply(model_obj$ranks[[rank]]$presets, function(pr) as.numeric(pr$x)))
  preset_df <- as.data.frame(preset_mat)
  names(preset_df) <- p$vnames
  cases_df <- rbind(p$data[, p$vnames, drop = FALSE], preset_df, rand)
  cases_json <- toJSON(lapply(seq_len(nrow(cases_df)), function(i) as.numeric(cases_df[i, ])), digits = 17, auto_unbox = TRUE)
  back <- fromJSON(cases_json)
  if (is.list(back) && !is.matrix(back)) back <- do.call(rbind, back)
  back <- as.matrix(back)
  if (!identical(as.numeric(back), as.numeric(as.matrix(cases_df)))) {
    stop(rank, ": inputs did not round-trip through JSON")
  }
  cases_df <- as.data.frame(back)
  names(cases_df) <- p$vnames
  rv <- predict(p$rf, cases_df, type = "vote", norm.votes = FALSE)
  rv <- unname(as.matrix(rv[, p$classes, drop = FALSE]))
  r_class <- as.character(predict(p$model, cases_df))
  rf_class <- as.character(predict(p$rf, cases_df))
  one <- list(
    classes = I(unname(p$classes)),
    cutoff = I(unname(as.numeric(p$rf$forest$cutoff))),
    features = model_obj$ranks[[rank]]$features,
    trees = p$trees,
    importance = model_obj$ranks[[rank]]$importance
  )
  one_json <- toJSON(one, auto_unbox = FALSE, digits = 17, null = "null")
  invisible(ctx$eval(paste0("var model = ", one_json, ";")))
  invisible(ctx$eval("hydrateModel(model);"))
  invisible(ctx$eval(paste0("var cases = ", cases_json, ";")))
  invisible(ctx$eval("
var voteRows = new Array(cases.length);
var i, j, voteSum = 0;
for (i = 0; i < cases.length; i++) {
  voteRows[i] = predictCase(model, cases[i]).votes;
  for (j = 0; j < voteRows[i].length; j++) voteSum += voteRows[i][j];
}
var scaled = scaleImportance(model.importance.rawValues);
var scaledDiff = 0;
for (i = 0; i < scaled.length; i++) {
  var d = Math.abs(scaled[i] - model.importance.scaled[i]);
  if (d > scaledDiff) scaledDiff = d;
}
var payload = JSON.stringify({votes: voteRows, voteSum: voteSum, scaledDiff: scaledDiff});
"))
  payload <- fromJSON(ctx$get("payload"))
  js_votes <- payload$votes
  if (is.list(js_votes) && !is.matrix(js_votes)) js_votes <- do.call(rbind, js_votes)
  js_votes <- as.matrix(js_votes)
  if (!identical(dim(js_votes), dim(rv))) stop(rank, ": vote matrix size mismatch")
  if (max(abs(js_votes - rv)) > 0) stop(rank, ": JavaScript votes differ from R")
  cutoff <- as.numeric(p$rf$forest$cutoff)
  equal_cutoff <- all(cutoff == cutoff[[1]])
  ties <- 0L
  non_tie_bad <- 0L
  tie_bad <- 0L
  for (i in seq_len(nrow(js_votes))) {
    leaders <- if (equal_cutoff) {
      p$classes[js_votes[i, ] == max(js_votes[i, ])]
    } else {
      score <- js_votes[i, ] / cutoff
      p$classes[score == max(score)]
    }
    if (length(leaders) > 1) {
      ties <- ties + 1L
      if (!(r_class[i] %in% leaders) || !(rf_class[i] %in% leaders)) tie_bad <- tie_bad + 1L
    } else if (leaders[[1]] != r_class[i] || leaders[[1]] != rf_class[i]) {
      non_tie_bad <- non_tie_bad + 1L
    }
  }
  if (non_tie_bad > 0 || tie_bad > 0) {
    stop(rank, ": class mismatch non-tie ", non_tie_bad, " tie ", tie_bad)
  }
  message(rank, ": ", nrow(cases_df), " rows, votes match, ties ", ties)
  verify_ranks[[rank]] <- list(
    rows = unbox(nrow(cases_df)),
    trainingRows = unbox(nrow(p$data)),
    presetRows = unbox(nrow(preset_df)),
    randomRows = unbox(nrow(rand)),
    voteMismatches = unbox(0L),
    nonTieClassMismatches = unbox(non_tie_bad),
    ties = unbox(ties),
    tieOutsideLeaders = unbox(tie_bad),
    scaledDiff = unbox(payload$scaledDiff),
    equalCutoff = unbox(equal_cutoff)
  )
}

verify <- list(ranks = verify_ranks)
model_json <- toJSON(model_obj, auto_unbox = FALSE, digits = 17, null = "null", na = "null")
verify_json <- toJSON(verify, auto_unbox = FALSE, digits = 17, null = "null", na = "null")
model_json <- gsub("<", "\\u003c", model_json, fixed = TRUE)
verify_json <- gsub("<", "\\u003c", verify_json, fixed = TRUE)

insert_token <- function(html, token, value) {
  parts <- strsplit(html, token, fixed = TRUE)[[1]]
  if (length(parts) != 2L) stop("Placeholder ", token, " count is ", length(parts) - 1L)
  paste0(parts[[1]], value, parts[[2]])
}
template <- paste(readLines(file.path(root, "index.template.html"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
html <- insert_token(template, "___PREDICT_JS___", predict_src)
html <- insert_token(html, "___MODEL_JSON___", model_json)
html <- insert_token(html, "___VERIFY_JSON___", verify_json)
out <- file(file.path(root, "index.html"), open = "wb")
writeBin(charToRaw(enc2utf8(html)), out)
close(out)
message("Wrote index.html")
