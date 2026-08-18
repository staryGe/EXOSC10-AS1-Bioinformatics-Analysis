############# 1. 环境准备与包加载 #############

# 确保 pak 安装好
if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

# 检查并安装缺失的包
packages <- c("readxl", "tidyverse", "DESeq2")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    pak::pkg_install(pkg)
  }
}

library(readxl)
library(tidyverse)
library(DESeq2)

############# 2. 数据读取与矩阵预处理 #############

# 使用相对路径读取原始 Counts 数据
file_path <- "./data/GSE289213_Genes_expression.xlsx"
raw_df <- read_excel(file_path)

# 查看数据维度与前几行
dim(raw_df)
head(raw_df)

# 校验：检查除第一列外的数值是否均为整数 counts
stopifnot(all(raw_df[, -1] %% 1 == 0))

# 转换为表达矩阵，并将 gene_id 设为行名
counts_matrix <- raw_df %>% 
  column_to_rownames(var = "gene_id") %>% 
  as.matrix()

# 强制转换数值类型为整数 (integer)
mode(counts_matrix) <- "integer"

############# 3. 构建 DESeq2 容器与过滤 #############

# 构建样本表型信息表 (colData)
sample_info <- data.frame(
  row.names = colnames(counts_matrix),
  group = factor(c(
    rep("NC", 3),
    rep("siRNA_1", 3),
    rep("siRNA_2", 3)
  ), levels = c("NC", "siRNA_1", "siRNA_2")) # NC 作为 reference 对照组
)

# 校验行名与列名匹配
stopifnot(all(rownames(sample_info) == colnames(counts_matrix)))

# 创建 DESeq2 核心数据对象
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = sample_info,
  design    = ~ group
)

# 过滤低表达基因（保留在至少 3 个样本中 count >= 10 的基因）
keep <- rowSums(counts(dds) >= 10) >= 3
dds_filtered <- dds[keep, ]

cat("过滤前基因数:", nrow(dds), "\n过滤后基因数:", nrow(dds_filtered), "\n")

############# 4. 质量控制 (QC) 与 PCA 绘图 #############

# VST 快速标准化
vsd <- vst(dds_filtered, blind = FALSE)

# 绘制 PCA 质检图并自动保存
p_pca <- plotPCA(vsd, intgroup = "group") + 
  theme_bw() + 
  labs(title = "PCA Plot: Sample Quality Control")

print(p_pca)

# 保存质检图片至 figures 文件夹
if(!dir.exists("./figures")) dir.create("./figures")
ggsave("./figures/001_QC_PCA.png", plot = p_pca, width = 6, height = 5, dpi = 300)