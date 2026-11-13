# 微生物菌相：phyloseq 拆分類階層，再用 Random Forest 做 10-fold CV
# 與 R_Machine_learning_IRIS/01_machine_learning.R 相同的 caret 設定：
#   set.seed(123)、trainControl(method = "cv", number = 10)、method = "rf"

source(file.path(script_dir_01 <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }
}), "microbiota_common.R"))
root <- script_dir_01

bundle <- build_models(root)

message("\nBaseline (always predict the larger group): ", bundle$baseline_label)
acc_df <- data.frame(
  Rank = names(bundle$accuracy),
  Accuracy = as.numeric(bundle$accuracy),
  stringsAsFactors = FALSE
)
acc_df$Rank <- factor(acc_df$Rank, levels = RANKS)
print(acc_df, digits = 4)
message("Highest CV accuracy: ", bundle$best_rank)

dir.create(file.path(root, "figures"), showWarnings = FALSE)

p_rank <- ggplot(acc_df, aes(x = Rank, y = Accuracy)) +
  geom_col(fill = "#1F6F6A", width = 0.7) +
  geom_hline(yintercept = bundle$baseline, linetype = "dashed", color = "#655d54") +
  geom_text(aes(label = sprintf("%.3f", Accuracy)), vjust = -0.4, size = 3.6) +
  ylim(0, 1) +
  labs(
    title = "10-fold CV accuracy by taxonomic rank",
    subtitle = "Dashed line: majority-class baseline",
    x = NULL,
    y = "Mean accuracy"
  ) +
  theme_minimal()
print(p_rank)
ggsave(file.path(root, "figures", "01_rank_accuracy.png"), p_rank, width = 7.2, height = 4.4, dpi = 140)

best <- bundle$ranks[[bundle$best_rank]]
folds <- best$resample
folds$Resample <- factor(folds$Resample, levels = unique(folds$Resample[order(folds$Resample)]))
p_fold <- ggplot(folds, aes(x = Resample, y = Accuracy)) +
  geom_col(fill = "#B2182B", width = 0.7) +
  geom_hline(yintercept = best$accuracy, linetype = "dashed", color = "#4a453f") +
  labs(
    title = paste("10-fold CV —", bundle$best_rank),
    subtitle = "Random forest accuracy per fold",
    x = "Fold",
    y = "Accuracy"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
print(p_fold)
ggsave(file.path(root, "figures", "01_best_rank_folds.png"), p_fold, width = 7.2, height = 4.4, dpi = 140)

message("Saved ml_models.rds, rank_features.csv, and figures/")
