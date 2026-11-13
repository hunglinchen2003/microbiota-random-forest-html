# 看每個分類階層的變數重要性與混淆矩陣。
# 模型來自 01_machine_learning.R 存下的 ml_models.rds；若還沒有，這裡會先訓練。

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
root <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(root, "microbiota_common.R"))

bundle <- load_models(root)
dir.create(file.path(root, "figures"), showWarnings = FALSE)

message("Majority-class baseline: ", bundle$baseline_label)
message("Best rank by CV accuracy: ", bundle$best_rank)

for (rank in RANKS) {
  obj <- bundle$ranks[[rank]]
  cat("\n========== ", rank, " (", RANK_ZH[[rank]], ") ==========\n", sep = "")
  cat("tax_glom taxa:", obj$n_glom,
      " | used in the model:", obj$n_features,
      " | prevalence cutoff:", obj$prevalence_cutoff,
      " | capped at", MAX_FEATURES, ":", obj$capped, "\n")
  cat("mtry:", obj$best_mtry, " ntree:", obj$ntree, "\n")
  cat("Mean CV accuracy:", obj$accuracy_label, "\n")
  cat("Features:\n")
  print(obj$feature_log[, c("taxon", "prevalence", "mean_percent")], digits = 3, row.names = FALSE)

  pred <- predict(obj$model, obj$data)
  cat("\nTraining-set confusion matrix (not the CV accuracy):\n")
  print(confusionMatrix(pred, obj$data$Group))

  imp_df <- data.frame(
    Taxon = unname(obj$labels[names(obj$importance_scaled)]),
    Importance = as.numeric(obj$importance_scaled),
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[order(imp_df$Importance), , drop = FALSE]
  p <- ggplot(imp_df, aes(x = reorder(Taxon, Importance), y = Importance)) +
    geom_col(fill = "#B2182B") +
    coord_flip() +
    labs(
      title = paste("Variable importance —", rank),
      x = NULL,
      y = "Scaled importance (0-100)"
    ) +
    theme_minimal()
  ggsave(
    file.path(root, "figures", paste0("02_importance_", rank, ".png")),
    p, width = 7.4, height = 4.8, dpi = 140
  )
  if (rank == bundle$best_rank) print(p)
}

message("\nImportance plots are in figures/")
