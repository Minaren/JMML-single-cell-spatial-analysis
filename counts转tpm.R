# ===============================
# 1. 安装和加载所需包
# ===============================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DGEobj.utils", ask = FALSE, update = TRUE)
library(DGEobj.utils)

# ===============================
# 2. 读取 count 矩阵
# ===============================
counts <- read.csv("counts_raw_data.csv", header = TRUE, check.names = FALSE)

# 处理重复基因名
counts[,1] <- make.unique(as.character(counts[,1]))
rownames(counts) <- counts[,1]
counts <- counts[,-1]  # 删除第一列，只保留数值列

# 转为数值矩阵
counts_numeric <- as.matrix(sapply(counts, as.numeric))
rownames(counts_numeric) <- rownames(counts)

# NA counts 替换为 0
counts_numeric[is.na(counts_numeric)] <- 0

# ===============================
# 3. 删除全零样本列
# ===============================
counts_numeric <- counts_numeric[, colSums(counts_numeric) > 0]

# ===============================
# 4. 设置基因长度
# ===============================
default_length <- 1000  # 可根据需要修改
gene_length_vec <- rep(default_length, nrow(counts_numeric))
names(gene_length_vec) <- rownames(counts_numeric)

# ===============================
# 5. 计算 TPM
# ===============================
tpm <- convertCounts(
  countsMatrix = counts_numeric,
  unit = "TPM",
  geneLength = gene_length_vec,
  log = FALSE,
  normalize = "none"
)

# 查看前 5 行
head(tpm, 5)

# ===============================
# 计算 FPKM
# ===============================
fpkm <- convertCounts(
  countsMatrix = counts_numeric,
  unit = "FPKM",          # 这里改成 FPKM
  geneLength = gene_length_vec,
  log = FALSE,            # 设置 TRUE 可返回 log2(FPKM)
  normalize = "none"      # 可选 "TMM" 等归一化方法
)

# 查看前 5 行
head(fpkm, 5)

# 6. 保存结果
# ===============================
write.table(tpm, file = "TPM_matrix.txt", sep = "\t", quote = FALSE)

# 保存结果
write.csv(tpm, file = "TPM_matrix.csv", quote = FALSE)
write.csv(fpkm, file = "FPKM_matrix.csv", quote = FALSE)
