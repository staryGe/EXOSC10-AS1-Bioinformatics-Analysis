############# 1. 环境准备与包加载 #############

if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

packages <- c("clusterProfiler", "org.Hs.eg.db", "tidyverse")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    pak::pkg_install(pkg)
  }
}

# 使用 pak 专门检查并安装 pathview (Bioconductor 包)
if (!requireNamespace("pathview", quietly = TRUE)) {
  pak::pkg_install("bioc::pathview")
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(pathview)

############# 2. 读取差异表达数据并转换 ID #############

if (!exists("res_df_si1")) {
  if (file.exists("./data/DEG_results_siRNA1_vs_NC.csv")) {
    res_df_si1 <- read.csv("./data/DEG_results_siRNA1_vs_NC.csv")
  } else {
    stop("未找到 DEG 结果数据，请先运行 02_volcano_and_top10.R 脚本！")
  }
}

# 提取下调基因 ENSEMBL ID
down_genes <- res_df_si1 %>% filter(change == "Down") %>% pull(gene_id)

# 将 ENSEMBL ID 转换为 KEGG 认可的 ENTREZID
down_entrez <- bitr(
  geneID   = down_genes, 
  fromType = "ENSEMBL", 
  toType   = "ENTREZID", 
  OrgDb    = org.Hs.eg.db
)$ENTREZID

cat("成功转换为 ENTREZID 的下调基因数:", length(down_entrez), "\n")

############# 3. KEGG 富集分析与气泡图保存 #############

kegg_down <- enrichKEGG(
  gene          = down_entrez,
  organism      = "hsa",           # hsa 代表 Homo sapiens (人类)
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# 绘制气泡图
p_kegg_down <- dotplot(kegg_down, showCategory = 10) + 
  theme_bw() + 
  labs(
    title = "KEGG Pathways: Down-regulated Genes", 
    x = "Gene Ratio"
  )

print(p_kegg_down)

# 自动保存 KEGG 气泡图
if(!dir.exists("./figures")) dir.create("./figures")
ggsave("./figures/007_KEGG_down.png", plot = p_kegg_down, width = 6, height = 5, dpi = 300)

# 保存 KEGG 富集结果数据表
if(!dir.exists("./data")) dir.create("./data")
write.csv(as.data.frame(kegg_down), "./data/KEGG_down_results.csv", row.names = FALSE)

############# 4. pathview 官方信号通路图高亮绘制 #############

# 构造表达量向量（下调赋值 -1）
foldchanges <- rep(-1, length(down_entrez))
names(foldchanges) <- down_entrez

# 运行 pathview 绘制 Cytokine-cytokine receptor interaction (hsa04060)
pv.out <- pathview(
  gene.data  = foldchanges, 
  pathway.id = "hsa04060",
  species    = "hsa",
  limit      = list(gene = c(-2, 2)), 
  bins       = list(gene = 10)
)

# 将生成的 pathview 图片移动并重命名至 figures/ 目录中
if (file.exists("hsa04060.pathview.png")) {
  file.rename("hsa04060.pathview.png", "./figures/008-1_hsa04060.pathview.png")
}
if (file.exists("hsa04060.png")) {
  file.rename("hsa04060.png", "./figures/008-2_hsa04060.png")
}

cat("脚本 04 运行完成！KEGG 气泡图及 pathview 标记通路图已成功整理至 figures/ 文件夹。\n")