rm(list=ls())
getwd()
setwd("E:/")

# 读取 CSV 文件
df <- read.csv("growth_curve.csv")

# 确保数据格式正确
df$Time <- as.numeric(df$Time)
df$Strain <- as.factor(df$Strain)
df$Replicate <- as.factor(df$Replicate)

# 查看数据结构
head(df)

# 载入 ggplot2
library(ggplot2)
library(dplyr)

# 计算每个时间点、每个菌株的均值和标准差
df_summary <- df %>%
  group_by(Time, Strain) %>%
  summarise(Mean_OD = mean(OD600), SD_OD = sd(OD600), .groups = "drop")

# 创建菌株名称映射（使用斜体）
strain_labels <- c(
  "wt" = "italic('Wild type')",
  "n" = "italic('ΔnfxB')",
  "nc" = "italic('ΔnfxBΔmexC')",
  "nd" = "italic('ΔnfxBΔmexD')",
  "nj" = "italic('ΔnfxBΔoprJ')"
)




ggplot(df_summary, aes(x = Time, y = Mean_OD, group = Strain, color = Strain, fill = Strain)) +
  geom_line(size = 1) +  # 平均生长曲线
  geom_point(size = 2) +   # 数据点
  geom_errorbar(aes(ymin = Mean_OD - SD_OD, ymax = Mean_OD + SD_OD), width = 0.2, size = 0.8) +  # 添加误差棒
  #geom_ribbon(aes(ymin = Mean_OD - SD_OD, ymax = Mean_OD + SD_OD), alpha = 0.05) +  # 误差带
  labs(x = "Time(h)", y = "OD600") +
  scale_color_manual(name = "Strains", values = c("#FF0000", "#00FF00", "#000000", "#FF00FF", "#A020F0"),
                     labels = lapply(strain_labels, function(x) parse(text = x))) +
  scale_fill_manual(name = "Strains", values = c("#FF0000", "#00FF00", "#000000", "#FF00FF", "#A020F0"),
                    labels = lapply(strain_labels, function(x) parse(text = x))) +
  guides(fill = "none") +  # **去掉右下角的 fill 图例**
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    axis.line = element_line(color = "black")   # 添加坐标轴线
   
  )




