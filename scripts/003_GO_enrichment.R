############# 1. 环境准备与包加载 #############

if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

packages <- c("clusterProfiler", "org.Hs.eg.db", "tidyverse")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    pak::pkg_install(pkg)
  }
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

############# 2. 读取差异表达基因数据 #############

if (!exists("res_df_si1")) {
  if (file.exists("./data/DEG_results_siRNA1_vs_NC.csv")) {
    res_df_si1 <- read.csv("./data/DEG_results_siRNA1_vs_NC.csv")
  } else {
    stop("未找到 DEG 结果数据，请先运行 02_volcano_and_top10.R 脚本！")
  }
}

# 提取基因集
down_genes <- res_df_si1 %>% filter(change == "Down") %>% pull(gene_id)
up_genes   <- res_df_si1 %>% filter(change == "Up")   %>% pull(gene_id)
deg_genes  <- res_df_si1 %>% filter(change %in% c("Up", "Down")) %>% pull(gene_id)

cat("用于 GO 分析的基因数量 - 下调:", length(down_genes), "| 上调:", length(up_genes), "\n")

# 确保图片保存目录存在
if(!dir.exists("./figures")) dir.create("./figures")
if(!dir.exists("./data")) dir.create("./data")

############# 3. GO-BP 分析 (生物学过程) #############

# --- A. 纯下调基因 GO-BP ---
ego_bp_down <- enrichGO(
  gene          = down_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

p_bp_down <- dotplot(ego_bp_down, showCategory = 10) + 
  theme_bw() + 
  labs(title = "GO-BP: Down-regulated Genes (Suppressed Processes)", x = "Gene Ratio")

ggsave("./figures/004-1_GO_BP_down.png", plot = p_bp_down, width = 7, height = 6, dpi = 300)

# --- B. 全局/上调基因 GO-BP ---
ego_bp_all <- enrichGO(
  gene          = deg_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

p_bp_all <- dotplot(ego_bp_all, showCategory = 10) + 
  theme_bw() + 
  labs(title = "GO-BP: All DEGs Enrichment Overview", x = "Gene Ratio")

ggsave("./figures/004-2_GO_BP_all.png", plot = p_bp_all, width = 7, height = 6, dpi = 300)

############# 4. GO-CC 分析 (细胞成分 - 纯下调) #############

ego_cc_down <- enrichGO(
  gene          = down_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

p_cc_down <- dotplot(ego_cc_down, showCategory = 10) + 
  theme_bw() + 
  labs(title = "GO-CC: Down-regulated Genes (Location/Structure)", x = "Gene Ratio")

ggsave("./figures/005_GO_CC_down.png", plot = p_cc_down, width = 7, height = 6, dpi = 300)

############# 5. GO-MF 分析 (分子功能 - 纯下调) #############

ego_mf_down <- enrichGO(
  gene          = down_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

p_mf_down <- dotplot(ego_mf_down, showCategory = 10) + 
  theme_bw() + 
  labs(title = "GO-MF: Down-regulated Genes (Molecular Activity)", x = "Gene Ratio")

ggsave("./figures/006_GO_MF_down.png", plot = p_mf_down, width = 7, height = 6, dpi = 300)

############# 6. 导出 GO 富集结果数据表 #############

write.csv(as.data.frame(ego_bp_down), "./data/GO_BP_down_results.csv", row.names = FALSE)
write.csv(as.data.frame(ego_cc_down), "./data/GO_CC_down_results.csv", row.names = FALSE)
write.csv(as.data.frame(ego_mf_down), "./data/GO_MF_down_results.csv", row.names = FALSE)

cat("脚本 03 运行完成！所有 GO 富集气泡图及结果表已成功保存至 figures/ 与 data/ 文件夹。\n")