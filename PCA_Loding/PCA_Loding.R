#===============PCA/PCA 载荷 loading=============
#由于小样本数据利用随机森林得到的微生物重要性排序结果不理想
#OOB‑error = **1.0**，混淆矩阵全部预测错误；模型**完全没有区分对照组 / 水凝胶组的能力**
#得到每个门类：`MeanDecreaseAccuracy`、`MeanDecreaseGini`、`Permutation_P_value`（1000 次置换检验 P 值）
#注意看输出：绝大多数 Permutation_P_value >>0.05，仅有极少数接近 0.1，**没有任何类群达到传统显著性 α=0.05**。完全符合高维极小样本的表现。
#因此选择了PCA/PCA 载荷 loading作为小样本排序方法


library(readr)
library(dplyr)
abund_raw <- read_csv("HC.csv")
abund_raw_tbl <- as.data.frame(abund_raw)
rownames(abund_raw_tbl) <- abund_raw_tbl[["phylum"]]
abund_raw_tbl <- abund_raw_tbl[, !colnames(abund_raw_tbl) %in% "phylum"]
abund_raw_tbl[is.na(abund_raw_tbl)] <- 0
# 转置：行=样本，列=门类
abund_t <- t(abund_raw_tbl)
# PCA，注意：绝对丰度数值很大，建议做CLR变换（微生物群落常用，消除成分效应）
library(compositions)
abund_clr <- clr(abund_t + 1) # +1避免0值报错
pca_res <- prcomp(abund_clr, scale. = FALSE)
# 提取loading（门类对PCA轴的贡献，绝对值越高，越关键）
loading <- pca_res$rotation
loading_df <- as.data.frame(loading)
loading_df$phylum <- rownames(loading_df)
loading_df$PC1_abs <- abs(loading_df$PC1)
loading_df$PC2_abs <- abs(loading_df$PC2)
loading_df$sum_abs <- loading_df$PC1_abs + loading_df$PC2_abs
# 按PC1+PC2总载荷绝对值降序，得到群落内部驱动变异的top类群
loading_df <- loading_df %>% arrange(desc(sum_abs))
write_csv(loading_df,"HC_PCA_loading_keystone_taxa.csv")
print(head(loading_df,20))

#==================== 新增绘图代码 开始 ====================
library(ggplot2)
library(ggrepel)

# 提取PCA方差解释度，用于坐标轴标题
pca_var <- pca_res$sdev^2
pca_percent <- round(pca_var / sum(pca_var)*100, 1)

# --------图1：Top20门类 sum_abs 横向柱状图 色盲安全配色--------
top20_load <- head(loading_df,20)
# 因子排序，保证绘图顺序
top20_load$phylum <- factor(top20_load$phylum, levels = rev(top20_load$phylum))

# Okabe‑Ito 色盲安全色板，循环使用
okabe_ito <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000")

p1 <- ggplot(top20_load, aes(x = sum_abs, y = phylum, fill = phylum)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = rep(okabe_ito, 3)) +
  labs(x = "Sum of absolute loadings (|PC1|+|PC2|)",
       y = "Microbial phylum",
       title = "Top‑20 taxa contributing to HC community variation") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

ggsave("PCA_top20_loading_barplot.pdf", plot = p1, width = 8, height = 9)
print(p1)

# --------图2：PCA‑biplot 样本 + 门类载荷双标图 --------
# 提取样本坐标
sample_coord <- as.data.frame(pca_res$x[,c("PC1","PC2")])
sample_coord$SampleID <- rownames(sample_coord)

# 取top30门类用于biplot箭头，避免标签爆炸
top30_taxa <- head(loading_df,30)

p2 <- ggplot() +
  # 样本点（黑色，Okabe‑Ito黑色）
  geom_point(data = sample_coord, aes(x=PC1, y=PC2), size=3, color="#000000") +
  # 门类载荷箭头
  geom_segment(data = top30_taxa,
               aes(x=0,y=0,xend=PC1,yend=PC2),
               arrow = arrow(length = unit(0.2,"cm")),
               colour="#D55E00") +
  geom_text_repel(data = top30_taxa,
                  aes(x=PC1,y=PC2, label=phylum),
                  size = 2.7, max.overlaps = 40) +
  labs(x = paste0("PC1 (", pca_percent[1],"%)"),
       y = paste0("PC2 (", pca_percent[2],"%)"),
       title = "PCA biplot of HC microbial community (CLR‑transformed)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("PCA_biplot.pdf", plot = p2, width = 10, height = 8)
print(p2)
#==================== 新增绘图代码 结束 ====================

