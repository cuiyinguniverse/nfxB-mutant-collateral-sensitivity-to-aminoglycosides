#rm(list=ls())

install.packages("ggplot2")
library(ggplot2)
install.packages("dplyr")
library(dplyr)

install.packages("ggrepel")
library(ggrepel) # 用于防止标签重叠

install.packages("extrafont")
library(extrafont)
font_import()  # 导入系统字体
loadfonts(device = "win")  # 如果是 Windows 系统

install.packages("svglite")
library(svglite)



getwd()
setwd("E:/2023_12_20_mass_spectrometry_Proteomics_Databases_search_results/")


# 读取数据和感兴趣蛋白列表
data <- read.csv("./nfxB_wt_Difference.csv")
interest_proteins <- readLines("interest_proteins.txt")


# 筛选出感兴趣的蛋白数据
data <- data %>%
  mutate(Significance = case_when(
    log2FoldChange > 1 & Pvalue < 0.05 ~ "Up-regulated",
    log2FoldChange < -1 & Pvalue < 0.05 ~ "Down-regulated",
    TRUE ~ "Not-Significant"
  ))

# 标注感兴趣蛋白
data <- data %>%
  mutate(Is_Interest = ifelse(Protein_ID %in% interest_proteins, "Yes", "No"))


# 绘制火山图


p <- ggplot(data, aes(x = log2FoldChange, y = -log10(Pvalue), color = Significance)) +
  geom_point(alpha = 0.7) +   # 绘制所有点
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not-Significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +   # 添加竖线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # 添加水平线
  # 添加防止标签重叠的标签
  geom_label_repel(
    data = filter(data, Is_Interest == "Yes"),
    aes(label = Protein_ID),
    size = 3,  # 标签字体大小
    box.padding = 0.8,  # 标签与点的间距
    point.padding = 0.5,  # 点与引线的间距
    max.overlaps = 10,  # 限制最大重叠数
    segment.color = "black",  # 引线颜色
    label.size = 0.25,  # 标签边框的线条宽度
    fill = "white",  # 标签背景颜色
    color = "black"  # 标签文字颜色
  ) +
  labs(
    title = expression(italic(Delta) * italic("nfxB") ~ "vs wild type"),
    x = expression(Log[2] * "(Fold Change)"),   # 横坐标标签
    y = expression("-Log"[10] * "(P-Value)")    # 纵坐标标签
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),              # 去除背景网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1), # 添加边框
    legend.position = "right",                 # 图例放在右侧
    legend.title = element_blank(),            # 去掉图例标题
    legend.text = element_text(size = 10),     # 调整图例字体大小
    #text = element_text(family = "Arial"),    # 设置字体
    plot.title = element_text(
      hjust = 0.5,                             # 图标题居中
      size = 15,                               # 设置标题字体大小
      face = "bold.italic",                    # 设置标题为粗体和斜体
      color = "black"                          # 设置标题颜色
    ),
    axis.title.x = element_text(
      size = 12,                               # 横坐标标题字体大小
      hjust = 0.5,                             # 横坐标标题居中
      face = "bold"                            # 横坐标标题为粗体
    ),
    axis.title.y = element_text(
      size = 12,                               # 纵坐标标题字体大小
      hjust = 0.5,                             # 纵坐标标题居中
      face = "bold"                            # 纵坐标标题为粗体
    ),
    axis.text = element_text(size = 10),       # 调整坐标轴刻度字体大小
    axis.ticks = element_line(size = 0.5),     # 设置坐标轴刻度线
    legend.key = element_blank()               # 去掉图例背景框
  )


# 确保 ggplot 图背景为纯白
p <- p + 
  theme(
    panel.background = element_rect(fill = "white", color = NA), # 图表背景为白色
    plot.background = element_rect(fill = "white", color = NA)  # 整体背景为白色
  )

# 保存为 tiff 格式，背景为纯白
ggsave("./delta_nfxB_vs_wt_volcano_plot.tiff", p, width = 7, height = 4.5, dpi = 300, bg = "white")
ggsave("./delta_nfxB_vs_wt_volcano_plot.svg", p, width = 7, height = 4.5)
ggsave("./delta_nfxB_vs_wt_volcano_plot.pdf", p, width = 7, height = 4.5)




