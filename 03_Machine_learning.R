# Shiny：選擇分類階層，調整相對豐度，預測 Control 或 Parkinson
# 模型與 01_machine_learning.R 相同（ml_models.rds）。

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
root <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(root, "microbiota_common.R"))
library(shiny)

bundle <- load_models(root)
specs <- lapply(bundle$ranks, slider_spec)
rank_choices <- stats::setNames(RANKS, paste0(RANKS, "（", RANK_ZH[RANKS], "）"))

acc_df <- data.frame(
  Rank = factor(names(bundle$accuracy), levels = RANKS),
  Accuracy = as.numeric(bundle$accuracy)
)

ui <- fluidPage(
  titlePanel("菌相階層預測 — Random Forest（Control / Parkinson）"),
  sidebarLayout(
    sidebarPanel(
      selectInput("rank", "分類階層", choices = rank_choices, selected = bundle$best_rank),
      uiOutput("slider_panel"),
      fluidRow(
        column(4, actionButton("reset_all", "全體中位數", style = "width:100%")),
        column(4, actionButton("reset_control", "Control", style = "width:100%")),
        column(4, actionButton("reset_pd", "Parkinson", style = "width:100%"))
      ),
      hr(),
      h4("Model summary"),
      htmlOutput("summary"),
      plotOutput("importance_plot", height = "260px")
    ),
    mainPanel(
      h3("預測組別"),
      verbatimTextOutput("prediction"),
      plotOutput("prob_plot", height = "240px"),
      hr(),
      h4("各階層 10-fold CV 準確率"),
      plotOutput("rank_plot", height = "260px"),
      hr(),
      h4("這個 app 在做什麼"),
      p("OTU 表先做成 phyloseq，再用 tax_glom 依門、綱、目、科、屬、種加總。每個 taxon 的數值是相對豐度（%），分母是該樣本在這個階層的全部讀數，不是只在滑桿上的那些 taxon。"),
      p("進模型前先拿掉名稱含 unidentified、uncultured、ambiguous、metagenome 的項目，並保留盛行率夠高、變異較大的 taxon。每個階層最多 12 個特徵，所以滑桿就是模型真正用到的全部變數。"),
      p("訓練設定與鳶尾花練習相同：set.seed(123)、10-fold CV、caret method = \"rf\"。虛線是永遠猜人數較多的那一組（majority baseline）。")
    )
  )
)

server <- function(input, output, session) {
  output$slider_panel <- renderUI({
    rank <- input$rank
    spec <- specs[[rank]]
    tagList(lapply(seq_along(spec), function(i) {
      f <- spec[[i]]
      sliderInput(
        paste0("s_", rank, "_", i),
        f$label,
        min = f$min, max = f$max, value = f$value, step = f$step
      )
    }))
  })

  set_medians <- function(group) {
    rank <- input$rank
    med <- if (is.null(group)) rank_medians(bundle$ranks[[rank]]) else rank_medians(bundle$ranks[[rank]], group)
    spec <- specs[[rank]]
    for (i in seq_along(spec)) {
      updateSliderInput(session, paste0("s_", rank, "_", i), value = unname(med[[spec[[i]]$column]]))
    }
  }
  observeEvent(input$reset_all, set_medians(NULL))
  observeEvent(input$reset_control, set_medians("Control"))
  observeEvent(input$reset_pd, set_medians("Parkinson"))

  current_row <- reactive({
    rank <- input$rank
    spec <- specs[[rank]]
    vals <- lapply(seq_along(spec), function(i) {
      value <- input[[paste0("s_", rank, "_", i)]]
      if (is.null(value)) spec[[i]]$value else value
    })
    newdata <- as.data.frame(vals, stringsAsFactors = FALSE)
    colnames(newdata) <- vapply(spec, function(f) f$column, character(1))
    newdata
  })

  pred_result <- reactive({
    newdata <- current_row()
    model <- bundle$ranks[[input$rank]]$model
    list(
      pred = predict(model, newdata),
      prob = predict(model, newdata, type = "prob")
    )
  })

  output$prediction <- renderText({
    res <- pred_result()
    prob <- res$prob[1, , drop = TRUE]
    paste0(
      "Predicted group: ", as.character(res$pred), "\n\n",
      paste(names(prob), round(as.numeric(prob), 3), collapse = "  |  ")
    )
  })

  output$prob_plot <- renderPlot({
    res <- pred_result()
    df <- data.frame(Group = names(res$prob), Probability = as.numeric(res$prob[1, ]))
    ggplot(df, aes(x = Group, y = Probability, fill = Group)) +
      geom_col(width = 0.65, show.legend = FALSE) +
      scale_fill_manual(values = c(Control = "#2A6F7F", Parkinson = "#C4622D")) +
      ylim(0, 1) +
      theme_minimal() +
      labs(title = paste("Prediction probability —", input$rank), x = NULL, y = "Probability")
  })

  output$summary <- renderUI({
    obj <- bundle$ranks[[input$rank]]
    HTML(paste0(
      "<p>Algorithm: Random Forest</p>",
      "<p>Validation: 10-fold Cross Validation</p>",
      "<p>Average Accuracy: ", obj$accuracy_label, "</p>",
      "<p>Majority baseline: ", bundle$baseline_label, "</p>",
      "<p>mtry: ", obj$best_mtry, " &nbsp; trees: ", obj$ntree, "</p>",
      "<p>tax_glom taxa: ", obj$n_glom, " &nbsp; in the model: ", obj$n_features, "</p>",
      "<p>Prevalence cutoff: ", obj$prevalence_cutoff, "</p>"
    ))
  })

  output$importance_plot <- renderPlot({
    obj <- bundle$ranks[[input$rank]]
    df <- data.frame(
      Taxon = unname(obj$labels[names(obj$importance_scaled)]),
      Importance = as.numeric(obj$importance_scaled)
    )
    ggplot(df, aes(x = reorder(Taxon, Importance), y = Importance)) +
      geom_col(fill = "#B2182B") +
      coord_flip() +
      theme_minimal() +
      labs(title = "Variable Importance (Random Forest)", x = NULL, y = "Importance")
  })

  output$rank_plot <- renderPlot({
    ggplot(acc_df, aes(x = Rank, y = Accuracy)) +
      geom_col(fill = "#1F6F6A", width = 0.7) +
      geom_hline(yintercept = bundle$baseline, linetype = "dashed") +
      ylim(0, 1) +
      theme_minimal() +
      labs(title = "CV accuracy by rank", subtitle = "Dashed line: majority baseline", x = NULL, y = "Mean accuracy")
  })
}

shinyApp(ui, server)
