server <- function(input, output, session) {
  # 观察 peri 事件-------------------------------------------
  observeEvent(input$action_peri, {
    # 载入此页面的建模数据
    dat <- readRDS("data/sarcoma_peri_train_dataset.rds")

    dd <<- datadist(dat)
    options(datadist = "dd")

    model <- cph(
      Surv(time, dead == 1) ~ age + grade2 + figo + chem + peri,
      x = T,
      y = T,
      data = dat,
      surv = T
    )
    nom <- nomogram(model, maxscale = 100)

    data_input <- data.frame(
      "age" = input$age_peri,
      "grade2" = input$grade2_peri,
      "figo" = input$figo_peri,
      "chem" = input$chem_peri,
      "peri" = input$peri
    )

    # 总结个案
    output$case_report_peri <- renderTable(
      {
        case <- data_input %>%
          mutate(
            grade2 = case_when(
              grade2 == 1 ~ "Low-grade",
              grade2 == 2 ~ "High-grade",
              grade2 == 3 ~ "Gx"
            ),
            figo = case_when(
              figo == 1 ~ "Stage I",
              figo == 2 ~ "Stage II",
              figo == 3 ~ "Stage III",
              figo == 4 ~ "Stage IV"
            ),
            chem = case_when(
              chem == 0 ~ "No",
              chem == 1 ~ "Yes"
            ),
            peri = case_when(
              peri == 0 ~ "Negtive",
              peri == 1 ~ "Malignant"
            )
          )
        colnames(case) <- c(
          "Age (years)", "Grade", "FIGO stage", "Adjuvant chemotherapy",
          "Peritoneal cytology status"
        )
        case
      },
      striped = T,
      hover = T,
      bordered = T,
      rownames = F,
      colnames = TRUE,
      align = "l",
      digits = 0
    )

    # 计算总分
    output$points_peri <- renderText({
      data_input[, 1:ncol(data_input)] <- lapply(data_input[, 1:ncol(data_input)], as.numeric)
      points <- points_cal(formula_rd(nom)[["formula"]], rd = data_input)
      print(round(points, digits = 1))
    })

    # 绘制预测生存曲线
    output$os_curve_peri <- renderPlot({
      surv_pred <- survest(
        model,
        newdata = data_input,
        conf.int = 0.95
      )

      # 转换为数据框以适配 ggplot2
      surv_df <- data.frame(
        time = as.numeric(surv_pred$time),
        survival = as.numeric(surv_pred$surv),
        lower = as.numeric(surv_pred$lower),
        upper = as.numeric(surv_pred$upper)
      )


      # 使用 ggplot2 进行绘制
      ggplot(surv_df, aes(x = time, y = survival)) +
        geom_line(linewidth = 1, color = "#E69A9C") + # 生存曲线
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, fill = "#E69A9C") + # 置信区间
        labs(
          x = "Time in months",
          y = "Survival probability (%)"
        ) +
        theme_minimal()
    })

    # 预测 OS
    output$predictions_peri <- renderTable(
      {
        survprob <- predictSurvProb(
          model,
          newd = data_input,
          times = c(1 * 12, 3 * 12, 5 * 12)
        )

        os_predict <- as.data.frame(t(survprob) * 100)

        row.names(os_predict) <- c(
          "Predicted 1-year OS (%)",
          "Predicted 3-year OS (%)",
          "Predicted 5-year OS (%)"
        )


        print(os_predict)
      },
      striped = T,
      hover = T,
      bordered = T,
      rownames = T,
      colnames = FALSE,
      align = "l",
      digits = 2
    )
  })









  # 观察 rad 事件-------------------------------------------------------------
  observeEvent(input$action_rad, {
    model <- readRDS("data/sarcoma_rad_model.rds")
    nom <- readRDS("data/sarcoma_rad_nomogram.rds")

    data_input <- data.frame(
      "age" = if_else(
        input$age < 49, 1,
        if_else(input$age <= 58, 2, 3)
      ) %>% as.factor(),
      "grade" = input$grade %>% as.factor(),
      "tumor_size" = input$tumor_size %>% as.factor(),
      "his" = input$his %>% as.factor(),
      "T_stage" = input$T_stage %>% as.factor(),
      "N_stage" = input$N_stage %>% as.factor(),
      "chemotherapy" = input$chemotherapy %>% as.factor()
    )

    points <- points_cal(formula_rd(nom)[["formula"]], rd = data_input)


    # 总结个案
    output$case_report_rad <- renderTable(
      {
        case <- data_input %>%
          mutate(
            age = case_when(
              age == 1 ~ "<49",
              age == 2 ~ "49–58",
              age == 3 ~ ">58"
            ),
            grade = case_when(
              grade == 1 ~ "Low-grade",
              grade == 2 ~ "High-grade",
              grade == 3 ~ "Gx"
            ),
            tumor_size = case_when(
              tumor_size == 1 ~ "<70 mm",
              tumor_size == 2 ~ "70–165 mm",
              tumor_size == 3 ~ ">165 mm",
              tumor_size == 4 ~ "Unspecific"
            ),
            his = case_when(
              his == 1 ~ "Leiomyosarcoma",
              his == 2 ~ "Endometrial stromal sarcoma",
              his == 3 ~ "Adenosarcoma"
            ),
            T_stage = case_when(
              T_stage == 1 ~ "T1",
              T_stage == 2 ~ "T2",
              T_stage == 3 ~ "T3",
              T_stage == 4 ~ "T4"
            ),
            N_stage = case_when(
              N_stage == 0 ~ "N0",
              N_stage == 1 ~ "N1"
            ),
            chemotherapy = case_when(
              chemotherapy == 0 ~ "No",
              chemotherapy == 1 ~ "Yes"
            )
          )
        colnames(case) <- c(
          "Age (years)", "Grade", "Tumor size", "Histology", "T stage",
          "N stage", "Adjuvant chemotherapy"
        )
        case
      },
      striped = T,
      hover = T,
      bordered = T,
      rownames = F,
      colnames = TRUE,
      digits = 0
    )

    output$points_rad <- renderText({
      print(round(points, digits = 1))
    })


    output$risk_rad <- renderUI({
      risk <- case_when(
        points < 135 ~ "Low-risk group",
        points < 221 ~ "Intermediate-risk group",
        points >= 221 ~ "High-risk group"
      )
      # 根据结果设置不同的背景颜色
      color <- case_when(
        points < 135 ~ "seagreen",
        points < 221 ~ "darkblue",
        points >= 221 ~ "#E69A9C"
      )

      HTML(
        paste(
          "<span style='color:",
          color,
          ";'>",
          risk,
          "</span>",
          sep = ""
        )
      )
    })

    output$benefit_rad <- renderUI({
      benefit <- if_else(
        points < 221, "Unable to achieve significant OS improvement",
        "Possible significant OS improvement"
      )
      # 根据结果设置不同的背景颜色
      color <- if_else(points < 221, "#E69A9C", "seagreen")

      HTML(
        paste(
          "<span style='color:",
          color,
          ";'>",
          benefit,
          "</span>",
          sep = ""
        )
      )
    })

    output$predictions_rad <- renderTable(
      {
        survprob <- predictSurvProb(
          model,
          newd = data_input,
          times = c(1 * 12, 3 * 12, 5 * 12)
        )

        os_predict <- as.data.frame(t(survprob) * 100)

        row.names(os_predict) <- c(
          "Predicted 1-year OS (%)",
          "Predicted 3-year OS (%)",
          "Predicted 5-year OS (%)"
        )


        print(os_predict)
      },
      striped = T,
      hover = T,
      bordered = T,
      rownames = T,
      colnames = FALSE,
      align = "l",
      digits = 2
    )
    output$os_curve_rad <- renderPlot({
      surv_pred <- survest(
        model,
        newdata = data_input,
        conf.int = 0.95
      )

      # 转换为数据框以适配 ggplot2
      surv_df <- data.frame(
        time = as.numeric(surv_pred$time),
        survival = as.numeric(surv_pred$surv),
        lower = as.numeric(surv_pred$lower),
        upper = as.numeric(surv_pred$upper)
      )


      # 使用 ggplot2 进行绘制
      ggplot(surv_df, aes(x = time, y = survival)) +
        geom_line(linewidth = 1, color = "#E69A9C") + # 生存曲线
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, fill = "#E69A9C") + # 置信区间
        labs(
          x = "Time in months",
          y = "Survival probability (%)"
        ) +
        theme_minimal()
    })
  })









  # 观察 ov事件-------------------------------------------------------------
  observeEvent(input$action_ov, {
    model <- readRDS("data/ov_sir_model.rds")
    nom <- readRDS("data/ov_nomogram.rds")

    data_input <- data.frame(
      "age" = input$age_ov %>% as.numeric(),
      "CA125" = input$CA125_ov %>% as.numeric(),
      "HE4" = input$HE4_ov %>% as.numeric(),
      "FIB" = input$FIB_ov %>% as.numeric(),
      "MLR" = if_else((input$mon_ov / input$lym_ov) > 0.25, 1, 0) %>% as.numeric(),
      "Bloodflow" = input$Bloodflow_ov %>% as.numeric(),
      "Solid_areas" = input$Solid_areas_ov %>% as.numeric()
    )

    data_case <- data.frame(
      "age" = input$age_ov,
      "CA125" = input$CA125_ov %>% factor(levels = c(0, 1), labels = c("≤ 100", "> 100")),
      "HE4" = input$HE4_ov %>% factor(levels = c(0, 1), labels = c("≤ 64.5", "> 64.5")),
      "FIB" = input$FIB_ov %>% factor(levels = c(0, 1), labels = c("≤ 3.28", "> 3.28")),
      "MLR" = if_else((input$mon_ov / input$lym_ov) > 0.25, 1, 0) %>% factor(levels = c(0, 1), labels = c("≤ 0.25", "> 0.25")),
      "Bloodflow" = input$Bloodflow_ov %>% factor(levels = c(0, 1), labels = c("No", "Yes")),
      "Solid_areas" = input$Solid_areas_ov %>% factor(levels = c(0, 1), labels = c("No", "Yes"))
    )

    malignant_prop <- predict(model, newdata = data_input, type = "response")
    malignant_prop <- round(malignant_prop * 100, 1)
    points <- points_cal(formula_rd(nom)[["formula"]], rd = data_input)


    # 总结个案
    output$case_report_ov <- renderTable(
      {
        data_case
      },
      striped = T,
      hover = T,
      bordered = T,
      rownames = F,
      colnames = TRUE,
      digits = 0
    )

    output$points_ov <- renderText({
      print(round(points, digits = 1))
    })


    output$predict_ov <- renderUI({
      # 根据结果设置不同的背景颜色
      color <- case_when(
        malignant_prop < 10 ~ "seagreen",
        malignant_prop < 30 ~ "darkblue",
        malignant_prop >= 30 ~ "#E69A9C"
      )

      HTML(
        paste(
          "<span style='color:",
          color,
          ";'>",
          malignant_prop,
          "%",
          "</span>",
          sep = ""
        )
      )
    })
  })









  # 观察 sc 事件--------------------------------------------------------
  observeEvent(input$action_sc, {
    # 载入此页面的建模数据
    gene_model <- readRDS("data/sc_gene_model.rds")
    final_model <- readRDS("data/sc_final_model.rds")



    gene_input <- data.frame(
      "C5AR1" = input$C5AR1_sc,
      "EREG" = input$EREG_sc,
      "FCGR2A" = input$FCGR2A_sc,
      "PLK3" = input$PLK3_sc,
      "SAT1" = input$SAT1_sc,
      "SLC11A2" = input$SLC11A2_sc,
      "SPAG9" = input$SPAG9_sc,
      "SRGAP1" = input$SRGAP1_sc,
      "TAOK3" = input$TAOK3_sc,
      "XBP1" = input$XBP1_sc
    )

    # 基因模型预测
    gene_score <- predict(gene_model, gene_input)$predicted

    clinic_input <- data.frame(
      "age_at_diagnosis" = input$age_sc,
      "residual_disease" = input$residual_sc %>%
        factor(
          levels = c(
            "No Macroscopic disease",
            "1-10 mm",
            "11-20 mm",
            ">20 mm"
          )
        ),
      "chemo_combined" = input$chemo_sc %>%
        factor(
          levels = c(
            "Paclitaxel + Carboplatin",
            "No chemotherapy",
            "Paclitaxel + Cisplatin",
            "More than 2 drugs",
            "Others"
          )
        ),
      "bevacizumab" = input$beva_sc %>% factor(levels = c("No", "Yes")),
      "gene_score" = gene_score
    )



    # 总结个案
    output$case_report_sc <- renderTable(
      {
        data.frame(
          "age_at_diagnosis" = input$age_sc,
          "residual_disease" = input$residual_sc,
          "chemo_combined" = input$chemo_sc,
          "bevacizumab" = input$beva_sc
        )
      },
      striped = T,
      hover = T,
      bordered = T,
      rownames = F,
      colnames = TRUE,
      align = "l",
      digits = 0
    )


    output$points_sc <- renderText({
      print(round(gene_score, digits = 1))
    })

    # 绘制预测生存曲线
    output$os_curve_sc <- renderPlot({
      surv_pred <- survest(
        final_model,
        newdata = clinic_input,
        conf.int = 0.95
      )

      # 转换为数据框以适配 ggplot2
      surv_df <- data.frame(
        time = as.numeric(surv_pred$time),
        survival = as.numeric(surv_pred$surv),
        lower = as.numeric(surv_pred$lower),
        upper = as.numeric(surv_pred$upper)
      )


      # 使用 ggplot2 进行绘制
      ggplot(surv_df, aes(x = time, y = survival)) +
        geom_line(linewidth = 1, color = "#E69A9C") + # 生存曲线
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, fill = "#E69A9C") + # 置信区间
        labs(
          x = "Time in months",
          y = "Survival probability (%)"
        ) +
        theme_minimal()
    })

    # 预测 OS
    output$predictions_sc <- renderTable(
      {
        survprob <- predictSurvProb(
          final_model,
          newd = clinic_input,
          times = c(1 * 12, 3 * 12, 5 * 12)
        )

        os_predict <- as.data.frame(t(survprob) * 100)

        row.names(os_predict) <- c(
          "Predicted 1-year OS (%)",
          "Predicted 3-year OS (%)",
          "Predicted 5-year OS (%)"
        )


        print(os_predict)
      },
      striped = T,
      hover = T,
      bordered = T,
      rownames = T,
      colnames = FALSE,
      align = "l",
      digits = 2
    )
  })
}
