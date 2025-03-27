{
  library(tidyverse)
  library(qs)
  library(Seurat)
  library(beepr)
  library(paletteer)
  library(cowplot)
  library(scales) # 用于自定义坐标轴刻度、显示颜色集
  library(decontX)
}


# 定义主题-------------------------------------------------------------------------------
# 条形图主题
theme_bar <- theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), # 去除主要和次要垂直网格线
    title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  )

# 密度图主题
theme_density <- theme_bw() +
  theme(
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# 小提琴图主题
theme_violin <- theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )


# 构建质量评价自编函数-------------------------------------------------------------------

sc.qc <- function(obj_name) {
  # 根据名称从环境中获取对应的Seurat对象
  seurat_obj <- get(obj_name)
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "【{obj_name}】信息\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()
  print(seurat_obj)

  # 判断是过滤前的质量评价还是过滤后（后缀带有“_filtered”）的质量评价
  is.filtered <- str_detect(obj_name, "filtered")

  if (is.filtered == FALSE) {
    # 建立文件夹，储存质控图像
    fig_path <<- str_glue("figures/QC/qc_{obj_name}/")
    if (!dir.exists(fig_path)) {
      dir.create(fig_path)
    }

    # 计算线粒体基因比例
    seurat_obj$mitoRatio <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    seurat_obj$mitoRatio <- seurat_obj$mitoRatio / 100

    str_glue(
      "---------------------------------------------------------------------------------\n",
      "【{obj_name}】的线粒体基因比例\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
    summary(seurat_obj$mitoRatio) %>% print()
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "质控指标已添加至【{obj_name}】的meta.data中\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
    head(seurat_obj) %>% print()

    # 将修改后的seurat_obj返回到全局环境中，覆盖原始的seurat_obj
    assign(obj_name, seurat_obj, envir = .GlobalEnv)
  }





  # 统计过滤前后各样本的细胞数用于绘图
  if (is.filtered) {
    # 从全局环境中获取过滤前的Seurat对象
    seurat_old <- get(str_remove(obj_name, "_filtered"))

    # 统计过滤前各样本的细胞数
    sample_old <- seurat_old@meta.data %>%
      group_by(sample_id) %>%
      summarise(count_old = n())

    # 统计过滤后各样本的细胞数
    sample_new <- seurat_obj@meta.data %>%
      group_by(sample_id) %>%
      summarise(count_new = n())

    # 合并过滤前后的样本数据
    QC_data <- left_join(
      sample_old,
      sample_new,
      by = "sample_id"
    ) %>%
      # 计算各个样本过滤掉了多少细胞
      mutate(
        diff = if_else(is.na(count_new), count_old, count_old - count_new)
      ) %>%
      # 转换成长数据
      pivot_longer(
        cols = c(count_new, diff),
        names_to = "cells",
        values_to = "count"
      ) %>%
      # 重命名标签
      mutate(
        cells = factor(
          cells,
          levels = c("diff", "count_new"),
          labels = c("Cells filtered", "Cells left")
        )
      )

    p_cells <- ggplot(QC_data, aes(x = sample_id, y = count, fill = cells)) +
      geom_bar(position = "stack", stat = "identity", width = .8) +
      scale_fill_manual(values = mycolors) +
      ggtitle(str_c("Number of cells of ", obj_name)) +
      xlab("Samples") +
      ylab("Number of cells") +
      scale_y_continuous(
        labels = label_comma(),
        expand = expansion(mult = c(0, 0.05)) # 让Y轴从0开始
      ) +
      theme_bar
    ggsave(width = 6, height = 4, str_glue("{fig_path}number_of_cells_after_filtering.pdf"))
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "细胞数条形图已保存为【{fig_path}number_of_cells_after_filtering.pdf】\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
  } else {
    if (length(unique(seurat_obj$sample_type)) > 1) {
      p_cells <- seurat_obj@meta.data %>%
        mutate(
          sample_id = factor(sample_id, levels = unique(sample_id[order(sample_type)]))
        ) %>%
        ggplot(aes(x = sample_id, fill = sample_type)) +
        geom_bar() +
        scale_fill_manual(values = mycolors) +
        ggtitle(str_c("Number of cells of ", obj_name)) +
        xlab("Samples") +
        ylab("Number of cells") +
        scale_y_continuous(
          labels = label_comma(),
          expand = expansion(mult = c(0, 0.05))
        ) +
        theme_bar
    } else {
      p_cells <- seurat_obj@meta.data %>%
        ggplot(aes(x = sample_id)) +
        geom_bar(fill = mycolors[1]) +
        ggtitle(str_c("Number of cells of ", obj_name)) +
        xlab("Samples") +
        ylab("Number of cells") +
        scale_y_continuous(
          labels = label_comma(),
          expand = expansion(mult = c(0, 0.05))
        ) +
        theme_bar
    }
    ggsave(width = 6, height = 4, str_glue("{fig_path}number_of_cells_original.pdf"))
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "细胞数条形图已保存为【{fig_path}number_of_cells_original.pdf】\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
  }
  print(p_cells)
  # 将p_cells重命名并返回到全局环境中
  assign(str_c(obj_name, "_p_cells"), p_cells, envir = .GlobalEnv)





  # Visualize the number UMIs(transcripts) per cell
  if(is.filtered) {
    p_umi <- seurat_obj@meta.data %>%
      ggplot(aes(x = nCount_RNA, color = sample_id, fill = sample_id)) +
      geom_density(alpha = 0.2) +
      scale_fill_manual(values = mycolors) +
      scale_color_manual(values = mycolors) +
      scale_x_log10(labels = label_comma()) + # 设置X轴刻度的小数格式
      geom_vline(xintercept = nCount_cut, lty = 2, color = "gray40") +
      annotate(
        "text",
        x = nCount_cut[1], y = 0, label = nCount_cut[1], # 标注阈值
        vjust = -0.5, hjust = -0.2, color = "black"
      ) +
      annotate(
        "text",
        x = nCount_cut[2], y = 0, label = nCount_cut[2], # 标注阈值
        vjust = -0.5, hjust = 1, color = "black"
      ) +
      xlab("Number of UMIs (log-scaled)") +
      ylab("Cell density") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      theme_density
  } else {
    p_umi <- seurat_obj@meta.data %>%
      ggplot(aes(x = nCount_RNA, color = sample_id, fill = sample_id)) +
      geom_density(alpha = 0.2) +
      scale_fill_manual(values = mycolors) +
      scale_color_manual(values = mycolors) +
      scale_x_continuous(labels = label_comma()) +
      geom_vline(xintercept = nCount_cut, lty = 2, color = "gray40") +
      xlab("Number of UMIs") +
      ylab("Cell density") +
      scale_y_continuous(
        labels = label_comma(),
        expand = expansion(mult = c(0, 0.05))
      ) +
      theme_density
  }



  # Visualize the distribution of genes detected per cell via histogram
  p_gene <- seurat_obj@meta.data %>%
    ggplot(aes(x = nFeature_RNA, color = sample_id, fill = sample_id)) +
    geom_density(alpha = 0.2) +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) +
    scale_x_continuous(labels = label_comma()) +
    geom_vline(xintercept = nFeature_cut, lty = 2, color = "gray40") +
    annotate(
      "text",
      x = nFeature_cut[1], y = 0, label = nFeature_cut[1], # 标注阈值
      vjust = -0.5, hjust = -0.2, color = "black"
    ) +
    annotate(
      "text",
      x = nFeature_cut[2], y = 0, label = nFeature_cut[2], # 标注阈值
      vjust = -0.5, hjust = 1.1, color = "black"
    ) +
    xlab("Number of genes") +
    ylab("Cell density") +
    scale_y_continuous(
      labels = label_comma(),
      expand = expansion(mult = c(0, 0.05))
      ) +
    theme_density

  # Visualize the distribution of mitochondrial gene expression detected per cell
  p_mito <- seurat_obj@meta.data %>%
    ggplot(aes(x = mitoRatio, color = sample_id, fill = sample_id)) +
    geom_density(alpha = 0.2) +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) +
    geom_vline(xintercept = mitoRatio_cut, lty = 2, color = "gray40") +
    scale_x_continuous(labels = label_number()) +
    annotate(
      "text",
      x = mitoRatio_cut, y = 0, label = str_c(mitoRatio_cut, "%"), # 标注阈值
      vjust = -0.5, hjust = 1.2, color = "black"
    ) +
    xlab("Proportion of mitochondrial genes (%)") +
    ylab("Cell density") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_density



  # 小提琴图
  violin_plot <- function(feature) {
    plot <- VlnPlot(
      seurat_obj,
      features = feature,
      group.by = "sample_id",
      split.by = "study_id",
      cols = mycolors,
      log = feature != "mitoRatio",
      pt.size = 0.1,
      alpha = 0.1
    ) +
      scale_y_continuous(labels = label_comma()) +
      geom_hline(
        yintercept = case_when(
          feature == "mitoRatio" ~ mitoRatio_cut,
          feature == "nCount_RNA" ~ nCount_cut,
          feature == "nFeature_RNA" ~ nFeature_cut
        ),
        lty = 2,
        color = "gray40"
      ) +
      xlab("Samples") +
      ylab(
        case_when(
          feature == "mitoRatio" ~ "Ratio of mitochondrial genes",
          feature == "nCount_RNA" ~ "Number of UMIs (log-scaled)",
          feature == "nFeature_RNA" ~ "Number of genes (log-scaled)"
        )
      ) +
      theme_violin
    return(plot)
  }


  p_vln <- plot_grid(
    violin_plot("nCount_RNA") + theme(legend.position = "none"),
    violin_plot("nFeature_RNA") + theme(legend.position = "none"),
    violin_plot("mitoRatio") + theme(legend.position = "none"),
    get_legend(violin_plot("mitoRatio")),
    nrow = 1,
    rel_widths = c(3, 3, 3, 1),
    labels = c("D", "E", "F")
  )



  p_qc <- plot_grid(
    p_umi + theme(legend.position = "none"),
    p_gene + theme(legend.position = "none"),
    p_mito + theme(legend.position = "none"),
    get_legend(p_mito),
    rel_widths = c(2, 2, 2, 1),
    nrow = 1,
    labels = "AUTO"
  ) %>%
    plot_grid(., p_vln, ncol = 1, rel_heights = c(1, 1.5))
  print(p_qc)

  if (is.filtered) {
    ggsave(width = 25, height = 11, str_glue("{fig_path}combined_QC_plot_after_filtering.pdf"))
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "合并质控图已保存为【{fig_path}combined_QC_plot_after_filtering.pdf】\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
  } else {
    ggsave(width = 25, height = 11, str_glue("{fig_path}combined_QC_plot_before_filtering.pdf"))
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "合并质控图已保存为【{fig_path}combined_QC_plot_before_filtering.pdf】\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
  }
  # 将p_qc重命名并返回到全局环境中
  assign(str_c(obj_name, "_p_qc"), p_qc, envir = .GlobalEnv)



  # 可视化线粒体基因含量和UMI的关系
  p_mito_umi <- seurat_obj@meta.data %>%
    ggplot(aes(x = nCount_RNA, y = mitoRatio, color = mitoRatio)) +
    geom_point(size = 0.1) +
    scale_color_gradient(low = mycolors[1], high = mycolors[2]) +
    geom_vline(xintercept = nCount_cut, lty = 2, color = "gray50") +
    geom_hline(yintercept = mitoRatio_cut, lty = 2, color = "gray50") +
    scale_x_continuous(labels = label_comma()) +
    scale_y_continuous(labels = label_comma()) +
    xlab("Number of UMIs") +
    ylab("Mitochondrial gene ratio") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    facet_wrap(~sample_id)
  print(p_mito_umi)

  if (is.filtered) {
    ggsave(
      width = 10.15,
      height = (length(unique(seurat_obj$sample_id)) / 4) %>% ceiling() * 2.45,
      str_glue("{fig_path}mito_vs_UMI_plot_after_filtering.pdf")
    )
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "线粒体基因含量和UMI的关系散点图已保存为【{fig_path}mito_vs_UMI_plot_after_filtering.pdf】\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
  } else {
    ggsave(
      width = 10.15,
      height = (length(unique(seurat_obj$sample_id)) / 4) %>% ceiling() * 2.45,
      str_glue("{fig_path}mito_vs_UMI_plot_before_filtering.pdf")
    )
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "线粒体基因含量和UMI的关系散点图已保存为【{fig_path}mito_vs_UMI_plot_before_filtering.pdf】\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
  }

  assign(str_c(obj_name, "_p_mito_umi"), p_mito_umi, envir = .GlobalEnv)


  # Joint filtering effects
  # Visualize the correlation between genes detected and number of UMIs
  # and determine whether strong presence of cells with low numbers of genes/UMIs
  p_join <- seurat_obj@meta.data %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) +
    geom_point(size = 0.2) +
    scale_colour_gradient(low = mycolors[1], high = mycolors[2]) +
    geom_smooth(method = lm, color = mycolors[3], lwd = 0.8) +
    scale_x_continuous(labels = label_comma()) +
    scale_y_continuous(labels = label_comma()) +
    geom_vline(xintercept = nCount_cut, lty = 2, color = "gray50") +
    geom_hline(yintercept = nFeature_cut, lty = 2, color = "gray50") +
    xlab("Number of UMIs") +
    ylab("Number of genes") +
    theme_bw() +
    facet_wrap(~sample_id)
  print(p_join)

  if (is.filtered) {
    ggsave(width = 12, str_glue("{fig_path}joint_QC_effects_after_filtering.pdf"))
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "联合质控指标散点图已保存为【{fig_path}joint_QC_effects_after_filtering.pdf】\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
  } else {
    ggsave(width = 12, str_glue("{fig_path}joint_QC_effects_before_filtering.pdf"))
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "联合质控指标散点图已保存为【{fig_path}joint_QC_effects_before_filtering.pdf】\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
  }

  assign(str_c(obj_name, "_p_join"), p_join, envir = .GlobalEnv)
}












# 构建细胞过滤自编函数--------------------------------------------------------------------

sc.filter <- function(obj_name, sample.filt = 0) {
  seurat_obj <- get(obj_name)
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "过滤前的【{obj_name}】\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()
  print(seurat_obj)
  seurat_filtered <- subset(
    seurat_obj,
    subset = (nCount_RNA >= nCount_cut[1]) &
      (nCount_RNA <= nCount_cut[2]) &
      (nFeature_RNA >= nFeature_cut[1]) &
      (nFeature_RNA <= nFeature_cut[2]) &
      (mitoRatio <= mitoRatio_cut) &
      !(seurat_obj$sample_id %in% sample.filt)
  )
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "质控过滤掉了【{ncol(seurat_obj) - ncol(seurat_filtered)}】个细胞\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "过滤后的【{obj_name}】\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()
  seurat_filtered %>% print()
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "过滤后各样本的细胞数量\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()
  table(seurat_filtered$sample_id) %>% print()
  assign(str_c(obj_name, "_filtered"), seurat_filtered, envir = .GlobalEnv)
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "过滤后的Seurat对象已保存至全局环境中，名称为【{obj_name}_filtered】\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()
}











# 构建评估和去除环境游离RNA污染自编函数---------------------------------------------------

sc.decont <- function(obj_name, decont.cut = 0.2, sample.filt = 0) {
  seurat_obj <- get(obj_name)
  print(seurat_obj)

  # 执行decontX，预测环境RNA比例
  if (str_detect(colnames(seurat_obj@meta.data), "Contamination") %>% any() == FALSE) {
    decont_result <- decontX(seurat_obj@assays$RNA$counts)
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "【{obj_name}】中的估计环境RNA污染情况\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
    summary(decont_result$contamination) %>% print()
    str_glue(
      "---------------------------------------------------------------------------------\n",
      "环境RNA百分比已添加至【{obj_name}】的mata.data中\n",
      "---------------------------------------------------------------------------------"
    ) %>% print()
    seurat_obj$Contamination <- decont_result$contamination
    assign(obj_name, seurat_obj, envir = .GlobalEnv)
    head(seurat_obj) %>% print()
  }


  p1 <- seurat_obj@meta.data %>%
    ggplot(aes(x = Contamination, color = sample_id, fill = sample_id)) +
    geom_density(alpha = 0.2) +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) +
    geom_vline(xintercept = decont.cut, lty = 2, color = "gray50") +
    annotate(
      "text",
      x = decont_cut, y = 0, label = decont_cut, # 标注阈值
      vjust = -0.5, hjust = -0.2, color = "black"
    ) +
    scale_x_continuous(labels = label_number()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    xlab("Contamination rate") +
    ylab("Cell density") +
    theme_density

  p2 <- VlnPlot(
    seurat_obj,
    features = "Contamination",
    group.by = "sample_id",
    split.by = "study_id",
    cols = mycolors,
    pt.size = 0.2,
    alpha = 0.2
  ) +
    geom_hline(yintercept = decont.cut, lty = 2, color = "gray40") +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    xlab("Samples") +
    ylab("Contamination rate") +
    theme_violin





  # 根据设定的阈值过滤高环境RNA污染的细胞
  seurat_low_con <- subset(
    seurat_obj,
    (Contamination <= decont.cut) &
      !(seurat_obj$sample_id %in% sample.filt)
  )
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "在设定阈值【{decont.cut}】下",
    "通过decontX过滤掉了【{ncol(seurat_obj) - ncol(seurat_low_con)}】个细胞\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()
  new_name <- str_split(obj_name, "_", simplify = TRUE)[, 1] %>%
    str_c(., "_low_con")
  assign(new_name, seurat_low_con, envir = .GlobalEnv)
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "过滤后的新的Seurat对象被命名为【{new_name}】并返回到全局环境中\n",
    "---------------------------------------------------------------------------------") %>%
    print()





  p3 <- seurat_low_con@meta.data %>%
    ggplot(aes(x = Contamination, color = sample_id, fill = sample_id)) +
    geom_density(alpha = 0.2) +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) +
    geom_vline(xintercept = decont.cut, lty = 2, color = "gray50") +
    annotate(
      "text",
      x = decont_cut, y = 0, label = decont_cut, # 标注阈值
      vjust = -0.5, hjust = -0.2, color = "black"
    ) +
    scale_x_continuous(labels = label_number()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    xlab("Contamination rate") +
    ylab("Cell density") +
    theme_density

  p4 <- VlnPlot(
    seurat_low_con,
    features = "Contamination",
    group.by = "sample_id",
    split.by = "study_id",
    cols = mycolors,
    pt.size = 0.2,
    alpha = 0.2
  ) +
    geom_hline(yintercept = decont.cut, lty = 2, color = "gray40") +
    scale_y_continuous(labels = label_number()) +
    xlab("Samples") +
    ylab("Contamination rate") +
    theme_violin

  p_decont <- plot_grid(
    p1, p2, p3, p4,
    ncol = 2,
    labels = "AUTO"
  )
  print(p_decont)
  ggsave(height = 10.1, width = 15.7, str_glue("{fig_path}decont.pdf"))
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "过滤前后的环境RNA分布图已保存至【{fig_path}decont.pdf】\n",
    "---------------------------------------------------------------------------------"
  ) %>%
    print()

  assign(str_c(obj_name, "_p_decont"), p_decont, envir = .GlobalEnv)






  # 绘制去除污染细胞前后各样本的细胞数统计条形图

  # 统计过滤前各样本的细胞数
  sample_old <- seurat_obj@meta.data %>%
    group_by(sample_id) %>%
    summarise(count_old = n())

  # 统计过滤后各样本的细胞数
  sample_new <- seurat_low_con@meta.data %>%
    group_by(sample_id) %>%
    summarise(count_new = n())

  # 合并过滤前后的样本数据
  data_plot <- left_join(
    sample_old,
    sample_new,
    by = "sample_id"
  ) %>%
    # 计算各个样本过滤掉了多少细胞
    mutate(
      diff = if_else(is.na(count_new), count_old, count_old - count_new)
    ) %>%
    # 转换成长数据
    pivot_longer(cols = c(count_new, diff), names_to = "cells", values_to = "count") %>%
    # 重命名标签
    mutate(
      cells = factor(
        cells,
        levels = c("diff", "count_new"),
        labels = c("Contaminated cells*", "Cells left")
      )
    )

  p_cells_after_decont <- ggplot(data_plot, aes(x = sample_id, y = count, fill = cells)) +
    geom_bar(position = "stack", stat = "identity", width = .8) +
    scale_fill_manual(values = mycolors) +
    ggtitle(str_glue("Number of cells after removal of contaminated cells in {obj_name}")) +
    xlab("Samples") +
    ylab("Number of cells") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # 让Y轴从0开始
    theme_bar
  print(p_cells_after_decont)
  ggsave(width = 6, height = 4, str_glue("{fig_path}number_of_cells_after_decont.pdf"))
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "去除污染细胞后的细胞数条形图已保存为【{fig_path}number_of_cells_after_decont.pdf】\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()

  assign(str_c(obj_name, "_p_cells_after_decont"), p_cells_after_decont, envir = .GlobalEnv)


  beep()
}


# 定义存储函数----------------------------------------------------------------------------

sc.save <- function(obj_name) {
  path <- str_split(obj_name, "_", simplify = TRUE)[, 1] %>%
    str_c("output/seurat_obj/", ., "_filtered.qs")
  qsave(get(obj_name), file = path)
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "过滤后的Seurat对象已保存为【{path}】\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()
}
