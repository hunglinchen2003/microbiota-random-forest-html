# 菌相階層的 Random Forest

OTU 表先用 phyloseq 彙整到門、綱、目、科、屬、種，再在每個階層分別做與鳶尾花練習相同的隨機森林。

網頁：<https://hunglinchen2003.github.io/microbiota-random-forest-html/>

直接打開 `index.html` 也可以離線使用。鳶尾花那一頁在 <https://hunglinchen2003.github.io/iris-random-forest-html/>。

## 資料

- `otu_table.csv`：OTU 讀數，菌在列、樣本在欄
- `taxonomy.csv`：Domain 到 Species
- `metadata.csv`：`Group` 為 Control 或 Parkinson

無法鑑定、ambiguous、uncultured、metagenome、unknown、incertae sedis 或空白的名稱，只在該階層改成 NA，再 `tax_glom(..., NArm = TRUE)`。相對豐度的分母是該階層的全部菌。盛行率依序嘗試 20%、10%、5%；若仍超過 12 個菌，留下相對豐度變異數最高的 12 個。拉桿使用的就是這一批菌。

## 執行

在這個資料夾、已安裝 phyloseq、caret、randomForest、ggplot2、shiny、jsonlite、V8 的 R 中：

```r
source("01_machine_learning.R")  # 訓練並存 ml_models.rds、圖
source("02_machine_learning.R")  # 各階層重要性與訓練集混淆矩陣
shiny::runApp("03_Machine_learning.R")
source("export_model.R")         # 產生 index.html
```

或用 Rscript 執行同名檔案。`01` 會在 `set.seed(123)` 之後立刻做 10-fold CV、`ntree = 500` 的 `caret::train(..., method = "rf")`。`index.html` 內嵌同一批樹；`predict.js` 的投票已在 R 的 V8 裡與 `predict(type = "vote", norm.votes = FALSE)` 核對。

平手時 R 會用亂數挑一個類別。網頁列出所有領先類別，不抽該亂數。
