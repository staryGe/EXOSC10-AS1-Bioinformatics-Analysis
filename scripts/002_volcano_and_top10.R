############# 1. 环境准备与包加载 #############

if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

packages <- c("tidyverse", "DESeq2", "org.Hs.eg.db", "pheatmap")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    pak::pkg_install(pkg)
  }
}

library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
library(pheatmap)

############# 2. 检查并获取上一步数据 #############

# 如果内存中没有 dds_filtered，尝试从本地 rds 加载
if (!exists("dds_filtered")) {
  if (file.exists("./data/dds_filtered.rds")) {
    dds_filtered <- readRDS("./data/dds_filtered.rds")
  } else {
    stop("未找到 dds_filtered 对象，请先运行 01_prepare_and_QC.R 脚本！")
  }
}

############# 3. 运行 DESeq2 核心计算 #############

# 运行 DESeq2 核心差异分析
dds <- DESeq(dds_filtered)

# 提取 siRNA_1 vs NC 的差异表达结果
res_siRNA1 <- results(dds, contrast = c("group", "siRNA_1", "NC"))

# 查看总结概要
summary(res_siRNA1)

############# 4. 差异基因筛选与注释 #############

# 转为标准 data.frame 并标记上/下调
res_df_si1 <- as.data.frame(res_siRNA1) %>% 
  rownames_to_column(var = "gene_id") %>% 
  mutate(
    change = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NOT"
    )
  )

# 输出差异基因统计
cat("--- 差异基因筛选统计 (padj < 0.05 & |log2FC| > 1) ---\n")
print(table(res_df_si1$change))

# 将 ENSEMBL ID 转换为 Gene Symbol
gene_map <- mapIds(
  org.Hs.eg.db,
  keys = res_df_si1$gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df_si1$symbol <- ifelse(is.na(gene_map), res_df_si1$gene_id, gene_map)

# 保存 DEG 筛选结果矩阵到 data/ 文件夹，供后续分析脚本调用
if(!dir.exists("./data")) dir.create("./data")
write.csv(res_df_si1, "./data/DEG_results_siRNA1_vs_NC.csv", row.names = FALSE)

############# 5. 绘制并保存火山图 #############

my_cols <- c("Down" = "#2f5597", "NOT" = "grey80", "Up" = "#c00000")

p_volcano <- ggplot(res_df_si1, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = my_cols) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_bw() +
  labs(
    title = "Volcano Plot: siRNA_1 vs NC",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-value)",
    color = "Significance"
  )

print(p_volcano)

# 自动导出火山图
if(!dir.exists("./figures")) dir.create("./figures")
ggsave("./figures/002_Volcano.png", plot = p_volcano, width = 6, height = 5, dpi = 300)

############# 6. 绘制并保存 Top10 基因热图 #############

# 提取上调 / 下调 Top10 基因
top10_up <- res_df_si1 %>% 
  filter(change == "Up") %>% 
  arrange(padj) %>% 
  head(10)

top10_down <- res_df_si1 %>% 
  filter(change == "Down") %>% 
  arrange(padj) %>% 
  head(10)

# 获取 VST 标准化表达矩阵
vsd <- vst(dds_filtered, blind = FALSE)
vsd_mat <- assay(vsd)

# --- A. 绘制并保存上调 Top 10 热图 ---
up_mat <- vsd_mat[top10_up$gene_id, ]
rownames(up_mat) <- top10_up$symbol

pheatmap(
  up_mat,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  main = "Top 10 Up-regulated Genes (siRNA_1 vs NC)",
  color = colorRampPalette(c("#2f5597", "white", "#c00000"))(100),
  filename = "./figures/003-1_up_Top10_heatmap.png", # 自动保存到图片文件夹
  width = 6,
  height = 5
)

# --- B. 绘制并保存下调 Top 10 热图 ---
down_mat <- vsd_mat[top10_down$gene_id, ]
rownames(down_mat) <- top10_down$symbol

pheatmap(
  down_mat,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  main = "Top 10 Down-regulated Genes (siRNA_1 vs NC)",
  color = colorRampPalette(c("#2f5597", "white", "#c00000"))(100),
  filename = "./figures/003-2_down_Top10_heatmap.png", # 自动保存到图片文件夹
  width = 6,
  height = 5
)

cat("脚本 02 运行完成！图表及 DEG 结果已成功保存。\n")