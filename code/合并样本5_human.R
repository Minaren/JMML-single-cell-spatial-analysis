############################################
# Figure 8 R 脚本（按图例顺序 A–K）
############################################

# --------------------------
# 准备工作
# --------------------------
rm(list = ls())
options(stringsAsFactors = F)
set.seed(220625)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))

# --------------------------
# 图 A：数据整合流程 & 质控起点
# --------------------------
dir = c('D:/aJMML/JMML/JMMLID5', 'D:/aJMML/JMML/PBM2')
samples_name= c("JMML","Healthy2")

scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                       min.cells=3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
  scRNAlist[[i]][["percent.MT"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
}

sc_merge <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])

# --------------------------
# 图 A：QC 可视化
# --------------------------
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.MT")
group = "orig.ident" 
VlnPlot(sc_merge, features = plot.featrures, group.by = group)

# QC 过滤
sc_filt <- subset(sc_merge, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.MT < 10)

# --------------------------
# 图 B：数据整合、UMAP 聚类
# --------------------------
sc_n = NormalizeData(sc_filt)
sc_n = FindVariableFeatures(sc_n, selection.method = "vst", nfeatures = 2000)
sc_n = ScaleData(sc_n)
sc_n = RunPCA(sc_n, npcs = 50)

sce_har = RunHarmony(sc_n, group.by.vars = "orig.ident")
sce.har = FindNeighbors(sce_har, dims = 1:30, reduction = "harmony")
sce.har = FindClusters(sce.har, resolution = 1.2, algorithm = 1)
sce.harm = RunUMAP(sce.har, dims = 1:30, reduction = "harmony")
DimPlot(sce.harm, label = T, repel = T)


# --------------------------
# 图 B：细胞类型注释 UMAP
# --------------------------
celltype = data.frame(ClusterID=0:19, celltype='NA')
celltype[celltype$ClusterID %in% c(9),2] = "HSC"
celltype[celltype$ClusterID %in% c(0),2] = "MPP"
celltype[celltype$ClusterID %in% c(1),2] = "CMP"
celltype[celltype$ClusterID %in% c(2),2] = "GMP"
celltype[celltype$ClusterID %in% c(3),2] = "MEP"
celltype[celltype$ClusterID %in% c(4),2] = "CLP"
celltype[celltype$ClusterID %in% c(5),2] = "B"
celltype[celltype$ClusterID %in% c(6),2] = "T"
celltype[celltype$ClusterID %in% c(7),2] = "NK"
celltype[celltype$ClusterID %in% c(8),2] = "Mono"
sce.harm$celltype = Idents(sce.harm)
DimPlot(sce.harm, label = T, repel = T)



# --------------------------
# 图 C：标记基因 DotPlot / Heatmap
# --------------------------
load("D:/aJMML/scRNA2/sce_har_anno.RData")
#1.读取目标genneset文件#
genesets <- read.csv("D:/aJMML/scRNA2/gsva_human_cluster.csv",header=F)
genesets <- subset(genesets, select = c("V1","V2")) %>% as.data.frame()
genesets <- split(genesets$V2, genesets$V1)

#2.提取分组平均表达矩阵#
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
expr <- AverageExpression(sce.harm, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res_mouse.csv", row.names = F)

my_colors <- colorRampPalette(c("#92b8da", "white", "#9f2b39"))(100) 
library(pheatmap)
# 绘制热图并应用自定义颜色
pheatmap(
  gsva.res,
  show_colnames = TRUE,
  scale = "row",  # 对行进行标准化
  cluster_row = FALSE,  # 不聚类行
  cluster_cols = FALSE,  # 不聚类列
  color = my_colors  # 应用自定义颜色
)



# --------------------------
# 图 E：细胞比例柱状图
# --------------------------
Cellratio <- prop.table(table(Idents(sce.harm), sce.harm$orig.ident), margin=2)
Cellratio_df <- as.data.frame(Cellratio)
ggplot(Cellratio_df, aes(x=Var2, y=Freq, fill=Var1)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(x="Sample", y="Cell proportion", fill="Cluster")

# --------------------------
# 图 D & E：HSC 亚群再分析
# --------------------------
HSC <- subset(sce.harm, idents = "HSC")
HSC <- ScaleData(HSC)
HSC <- RunPCA(HSC)
HSC <- FindNeighbors(HSC, dims=1:30)
HSC <- FindClusters(HSC, resolution=0.8)
HSC <- RunUMAP(HSC, dims=1:30)
DimPlot(HSC, label=T)

Cellratio_HSC <- prop.table(table(Idents(HSC), HSC$orig.ident), margin=2)
Cellratio_HSC_df <- as.data.frame(Cellratio_HSC)
ggplot(Cellratio_HSC_df, aes(x=Var2, y=Freq, fill=Var1)) +
  geom_bar(stat="identity") +
  theme_minimal()

# --------------------------
# 图 F：CD69 表达小提琴图 & FeaturePlot
# --------------------------
FeaturePlot(HSC, features="CD69", cols=c("lightgrey","#9f2b39"))
VlnPlot(HSC, features="CD69", pt.size=0, group.by="seurat_clusters") +
  geom_boxplot()

# --------------------------
# 图 G：细胞–细胞相互作用 & 增殖相关基因表达
# --------------------------
features_cell_interaction <- c("CD74","HLA-DRA","HLA-DPA1","HLA-DQB1","HLA-DQA1")
features_proliferation <- c("JUN","JUNB","FOS")
VlnPlot(HSC, features=features_cell_interaction)
VlnPlot(HSC, features=features_proliferation)

# --------------------------
# 图 H：增殖评分
# --------------------------
####增殖评分####
# 加载必要的包
library(Seurat)
library(ggplot2)
library(dplyr)

# 1. 人类增殖相关基因列表（对应你给的小鼠基因集改成了人类符号）
proliferation_gene_sets <- list(
  IEGs = c("FOS", "JUN", "JUNB", "JUND", "EGR1", "REL", "NFKBIA", "NFKBIE"),
  Cell_Cycle_Regulators = c("CCNG2", "CDK1", "CDK2", "CDK4", "CDK6", "CCNB1", "CCNB2",
                            "CCNA2", "CDC20", "CDC25A", "CDC25B", "CDC25C", "AURKB", "KIF23"),
  DNA_Replication_Mitosis = c("PCNA", "TOP2A", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
                              "UBE2C", "BIRC5", "TYMS", "RRM2", "CENPA", "CENPE", "CENPF"),
  Transcription_Factors = c("MYB", "FOXM1", "SOX4", "GATA2", "STAT5A", "ID1", "E2F1"),
  Signaling_Checkpoint_Others = c("TNF", "TNFAIP3", "JAK1", "KDM5B", "KDM6B",
                                  "CHEK1", "CHEK2", "CDK19")
)

# 合并所有增殖基因为一个向量，供AddModuleScore使用
proliferation_genes <- unique(unlist(proliferation_gene_sets))

# 2. 计算增殖分数
HSC <- AddModuleScore(
  object = HSC,
  features = list(proliferation_genes),
  name = "Proliferation_Score"
)

# 3. 可视化（小提琴图，默认按 seurat_clusters 分组）
VlnPlot(HSC, features = "Proliferation_Score1", group.by = "seurat_clusters", 
        pt.size = 0.1, cols = c("#9e2a2f","#4f85b8")) +
  labs(title = "Proliferation Score in HSC Clusters", y = "Proliferation Score") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("ProliferationScore_violin.png", width = 6, height = 5, dpi = 400)

# 4. Boxplot + 显著性检验示例
data <- FetchData(HSC, vars = c("Proliferation_Score1", "seurat_clusters"))
colnames(data) <- c("Score", "Cluster")
data$Cluster <- as.factor(data$Cluster)

ggplot(data, aes(x = Cluster, y = Score, fill = Cluster)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5) +
  labs(title = "Proliferation Score by HSC Cluster", x = "Cluster", y = "Proliferation Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("ProliferationScore_boxplot.png", width = 6, height = 5, dpi = 400)
ggsave("ProliferationScore_boxplot.pdf", width = 6, height = 5, dpi = 400)

# 5. 统计检验示例：cluster 0 vs cluster 1（你可根据实际cluster调整）
group0_vs_1 <- wilcox.test(Score ~ Cluster, data = subset(data, Cluster %in% c("0", "1")))
print("Wilcoxon test: Cluster 0 vs 1")
print(group0_vs_1)#p=1.215e-06

# --------------------------
# 图 I：GSVA 通路活性分析
# --------------------------
## 5. GSVA 分析
HSC_gsva <- gsva(expr = HSC_counts_filtered,
                 gset.idx.list = mouse_GO_bp_Set,
                 kcdf = "Poisson",
                 parallel.sz = 4)
write.table(HSC_gsva, "HSC_gsva_CD69.xls", sep = "\t")

## 6. 差异分析
group <- factor(meta[cells_to_use, "CD69_status"])
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(HSC_gsva, design)
cont.matrix <- makeContrasts(CD69high_vs_low = CD69high - CD69low, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, adjust = "fdr", number = Inf)
diff$ID <- rownames(diff)
write.csv(diff, "CD69_GSVA_diff.csv")

# 自定义感兴趣的通路
library(dplyr)
library(ggplot2)

# 确保diff中有ID列，如果没有，就添加
if (!"ID" %in% colnames(diff)) {
  diff$ID <- rownames(diff)
}

# 定义模块
proliferation_signaling_related <- c(
  "GOBP_POSITIVE_REGULATION_OF_ATP_BIOSYNTHETIC_PROCESS",
  "GOBP_POSITIVE_REGULATION_OF_NUCLEOTIDE_BIOSYNTHETIC_PROCESS",
  "GOBP_POSITIVE_REGULATION_OF_CELL_PROLIFERATION_IN_BONE_MARROW"
)

cell_interaction_related <- c(
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN",
  "GOBP_POSITIVE_REGULATION_OF_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_POSITIVE_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_CD4_POSITIVE_CD25_POSITIVE_ALPHA_BETA_REGULATORY_T_CELL_DIFFERENTIATION",
  "GOBP_REGULATION_OF_CD4_POSITIVE_CD25_POSITIVE_ALPHA_BETA_REGULATORY_T_CELL_DIFFERENTIATION"
)

up <- unique(c(proliferation_signaling_related, cell_interaction_related))

# 筛选并构造plot_df
plot_df <- diff %>%
  dplyr::filter(ID %in% up) %>%
  dplyr::select(ID, score = t) %>%  # 选择ID和t列，t重命名为score
  dplyr::mutate(module = dplyr::case_when(
    ID %in% proliferation_signaling_related ~ "Proliferation",
    ID %in% cell_interaction_related ~ "Cell Interaction",
    TRUE ~ "Other"
  )) %>%
  dplyr::arrange(module, score) %>%
  dplyr::mutate(ID = factor(ID, levels = unique(ID)))

# 画图
p <- ggplot(plot_df, aes(x = ID, y = score, fill = module)) +
  geom_bar(stat = 'identity', alpha = 0.9, color = "grey30", width = 0.7) +
  coord_flip() +
  facet_wrap(~ module, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold", color = "#4B4B4B"),
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "GSVA t value", title = "GSVA Pathways Grouped by Functional Module") +
  scale_fill_manual(values = c(
    "Proliferation" = "#9f2b39",
    "Cell Interaction" = "#4f85b8"
  ))

ggsave("gsva_barplot_faceted.pdf", plot = p, width = 12, height = 10)
# --------------------------
# 图 J：CD69high HSC-like 细胞示意（非代码生成）
# --------------------------
# 图 J 为示意图，显示 KRAS 突变 JMML 患者中 CD69high HSC-like 细胞丰度。

# --------------------------
# 图 K：Kaplan–Meier 生存分析
# ============================================================
# 完整 GSVA 方案：HSC_CD69high 活性评分 → 生存分析
# 一键运行（请先确保文件路径正确）
# ============================================================

# 加载包
library(Seurat)
library(GSVA)
library(dplyr)
library(survival)
library(survminer)
library(readxl)
library(data.table)

# ====================== 1. 数据准备 ======================

# --- 1.1 读取 bulk 表达矩阵 ---
# 文件名为 ids_exprs.csv，第一列为基因名
bulk_raw <- fread("ids_exprs.csv", data.table = FALSE)
rownames(bulk_raw) <- bulk_raw[, 1]          # 将第一列设为行名
bulk_raw <- bulk_raw[, -1, drop = FALSE]     # 去掉基因名列

# 转为数值矩阵，基因名大写，去重（取平均值）
bulk_expr_mat <- as.matrix(bulk_raw)
mode(bulk_expr_mat) <- "numeric"
rownames(bulk_expr_mat) <- toupper(rownames(bulk_expr_mat))

# 处理可能存在的重复基因
if (any(duplicated(rownames(bulk_expr_mat)))) {
  bulk_df <- as.data.frame(bulk_expr_mat)
  bulk_df$Gene <- rownames(bulk_df)
  bulk_df <- aggregate(. ~ Gene, data = bulk_df, FUN = mean)
  rownames(bulk_df) <- bulk_df$Gene
  bulk_df$Gene <- NULL
  bulk_expr_mat <- as.matrix(bulk_df)
  cat("Bulk 重复基因已去重，当前基因数：", nrow(bulk_expr_mat), "\n")
}

# --- 1.2 提取 HSC 亚群差异基因 ---
# 加载你保存的 HSC 对象（请确认路径）
load("D:/aJMML/scRNA2/human_HSC.RData")

# 设置亚群身份（0 = CD69high, 1 = CD69low）
Idents(HSC) <- "seurat_clusters"
HSC <- JoinLayers(HSC)   # Seurat v5 必须

# 寻找标记基因
hsc_markers <- FindAllMarkers(
  HSC,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)
write.csv(hsc_markers, "HSC_CD69high_vs_low_markers.csv", row.names = FALSE)
# 提取每个亚群的全部差异基因（已满足 logfc.threshold 和显著性）
genes_high <- hsc_markers %>%
  filter(cluster == "0") %>%
  pull(gene)

genes_low <- hsc_markers %>%
  filter(cluster == "1") %>%
  pull(gene)

# ====================== 2. GSVA 打分 ======================

common_high <- intersect(toupper(genes_high), rownames(bulk_expr_mat))
common_low  <- intersect(toupper(genes_low),  rownames(bulk_expr_mat))

cat("CD69high 基因交集数：", length(common_high), "\n")
cat("CD69low  基因交集数：", length(common_low), "\n")

gene_sets <- list(
  HSC_CD69high = common_high,
  HSC_CD69low  = common_low
)

# 兼容新旧 GSVA 版本
if (packageVersion("GSVA") >= "1.50.0") {
  param <- gsvaParam(as.matrix(bulk_expr_mat), gene_sets, kcdf = "Gaussian")
  gsva_res <- gsva(param, verbose = FALSE)
} else {
  gsva_res <- gsva(as.matrix(bulk_expr_mat), gene_sets, method = "ssgsea", kcdf = "Gaussian")
}

# 提取得分向量
hsc_high_score <- as.numeric(gsva_res["HSC_CD69high", ])
names(hsc_high_score) <- colnames(gsva_res)

# ====================== 3. 临床整合与生存分析 ======================

# --- 3.1 读取临床信息 ---
clinical <- read_excel("GSE71449_Table_S1.xlsx")   # 请根据实际文件名修改

# 统一 ID 列（假设原始 ID 列为 "ID"，否则修改）
clinical$SampleID <- as.character(clinical$ID)

# 处理生存变量（沿用你之前的字段名）
clinical$diag_survival_days <- as.numeric(clinical$`Survival from diagnosis\r\n(days)`)
clinical$diag_event <- ifelse(
  !is.na(clinical$`Cause of death`) & clinical$`Cause of death` != "", 1, 0
)

# --- 3.2 合并 GSVA 得分 ---
clinical$HSC_CD69high_score <- hsc_high_score[match(clinical$SampleID, names(hsc_high_score))]

# 去除缺失值
clinical_valid <- clinical[!is.na(clinical$HSC_CD69high_score) & 
                             !is.na(clinical$diag_survival_days) & 
                             !is.na(clinical$diag_event), ]

# --- 3.3 中位数分组 ---
med <- median(clinical_valid$HSC_CD69high_score, na.rm = TRUE)
clinical_valid$group_high <- ifelse(clinical_valid$HSC_CD69high_score > med, "High", "Low")
cat("中位数：", med, "\n分组情况：\n")
print(table(clinical_valid$group_high))

# --- 3.4 生存曲线 ---
surv_obj <- Surv(time = clinical_valid$diag_survival_days, event = clinical_valid$diag_event)
fit <- survfit(surv_obj ~ group_high, data = clinical_valid)
# 1. 绘制并赋值
km_plot <- ggsurvplot(fit, data = clinical_valid, pval = TRUE, risk.table = TRUE,
                      palette = "jco", xlab = "Days from diagnosis", ylab = "Overall Survival",
                      legend.title = "HSC_CD69high score")

# 2. 保存单独图形（不含风险表）
ggsave("KM_OS_HSC_CD69high.pdf", plot = km_plot$plot, width = 7, height = 6)

# 3. 保存完整图形（含风险表）
pdf("KM_OS_HSC_CD69high_full.pdf", width = 8, height = 7)
print(km_plot)
dev.off()



library(survival)
library(survminer)
# 1. 计算最佳截断点
cut <- surv_cutpoint(
  clinical_valid,
  time = "diag_survival_days",
  event = "diag_event",
  variables = "HSC_CD69high_score"
)

# 2. 根据截断值分组
clinical_valid$group_optimal <- ifelse(
  clinical_valid$HSC_CD69high_score > cut$cutpoint$cutpoint,
  "High", "Low"
)

# 3. 拟合生存曲线
surv_obj <- Surv(time = clinical_valid$diag_survival_days, event = clinical_valid$diag_event)
fit_opt <- survfit(surv_obj ~ group_optimal, data = clinical_valid)

# 4. 绘制 KM 曲线并保存
km_opt <- ggsurvplot(
  fit_opt, data = clinical_valid,
  pval = TRUE, risk.table = TRUE,
  palette = "jco", 
  xlab = "Days from diagnosis", ylab = "Overall Survival",
  legend.title = "HSC_CD69high (optimal cutoff)"
)

# 保存 PDF
pdf("KM_OS_HSC_CD69high_optimal_cutoff.pdf", width = 8, height = 7)
print(km_opt)
dev.off()

# 也可以只保存图形本身（不含风险表）：
# ggsave("KM_OS_HSC_CD69high_optimal_cutoff_plot.pdf", plot = km_opt$plot, width = 7, height = 6)

# ============================================================
# 完整流程：数据准备 → 单因素 → 多因素 → 森林图（方法一）
# ============================================================

# 加载包
library(survival)
library(survminer)
library(dplyr)
library(readxl)

# 1. 读取临床数据（请确认文件路径）
clinical <- read_excel("GSE71449_Table_S1.xlsx")   # 如文件名不同请修改
clinical$SampleID <- as.character(clinical$ID)      # 按实际 ID 列名调整

# 2. 处理生存终点
clinical$diag_survival_days <- as.numeric(clinical$`Survival from diagnosis\r\n(days)`)
clinical$diag_event <- ifelse(
  !is.na(clinical$`Cause of death`) & clinical$`Cause of death` != "", 1, 0
)

# 3. 读取 GSVA 得分（假设 hsc_high_score 已在工作环境中）
clinical$HSC_CD69high_score <- hsc_high_score[match(clinical$SampleID, names(hsc_high_score))]

# 4. 创建核心临床变量
clinical_valid <- clinical %>%
  filter(!is.na(HSC_CD69high_score) & !is.na(diag_survival_days) & !is.na(diag_event)) %>%
  mutate(
    age_years = as.numeric(`age at diagnosis\r\n(years)`),
    monosomy7 = ifelse(Karyotype == "mono7", 1, 0),
    mut_group = factor(case_when(
      Mutation == "PTPN11" ~ "PTPN11",
      Mutation == "quad neg" ~ "Quad_neg",
      TRUE ~ "Other"
    ), levels = c("PTPN11", "Other", "Quad_neg"))
  )

# 5. 单因素 Cox 回归（四个核心变量）
cat("========== 单因素分析 ==========\n")
vars <- c("HSC_CD69high_score", "age_years", "monosomy7", "mut_group")
for (v in vars) {
  f <- as.formula(paste("Surv(diag_survival_days, diag_event) ~", v))
  fit <- coxph(f, data = clinical_valid)
  cat("\n---", v, "---\n")
  print(summary(fit))
}
univ_table <- lapply(vars, function(v) {
  f <- as.formula(paste("Surv(diag_survival_days, diag_event) ~", v))
  fit <- coxph(f, data = clinical_valid)
  s <- summary(fit)
  coef <- s$coefficients
  ci <- s$conf.int
  
  # 数值变量：只有一行
  # 分类变量：每个水平一行
  data.frame(
    Variable = v,
    Level = rownames(coef),
    HR = round(coef[, "exp(coef)"], 2),
    Lower_95 = round(ci[, "lower .95"], 2),
    Upper_95 = round(ci[, "upper .95"], 2),
    P = signif(coef[, "Pr(>|z|)"], 3),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# 保存为 CSV
write.csv(univ_table, "Table_S_Univariate_Cox.csv", row.names = FALSE)
print(univ_table)



# 6. 多因素 Cox 回归（四个变量同时纳入）
cat("\n========== 多因素分析 ==========\n")
cox_multi <- coxph(
  Surv(diag_survival_days, diag_event) ~ HSC_CD69high_score + age_years + monosomy7 + mut_group,
  data = clinical_valid
)
summary(cox_multi)

# 7. 保存森林图（方法一：使用短变量名解决 ggforest 名称匹配错误）
# 1. 创建一个只包含模型所需短变量的纯数据框
df_plot <- data.frame(
  diag_survival_days = clinical_valid$diag_survival_days,
  diag_event         = clinical_valid$diag_event,
  score              = clinical_valid$HSC_CD69high_score,
  age                = clinical_valid$age_years,
  mono7              = clinical_valid$monosomy7,
  mut_group          = clinical_valid$mut_group   # 因子，水平为 PTPN11/Other/Quad_neg
)

# 彻底删除行名（可选）
rownames(df_plot) <- NULL

# 2. 用极简数据框重新拟合模型
cox_simple <- coxph(
  Surv(diag_survival_days, diag_event) ~ score + age + mono7 + mut_group,
  data = df_plot
)

# 3. 画森林图（这次数据框里只有这6列，匹配绝对不会出错）
library(survminer)
ggforest(cox_simple, data = df_plot)
pdf("Forest_plot_multivariate1.pdf", width = 9, height = 5)
ggforest(cox_simple, data = df_plot)
dev.off()



library(Seurat)
library(ggplot2)
library(patchwork)

# 你的两种颜色
group_colors <- c("#9e2a2f", "#4f85b8")

# 顶刊简约主题
top_theme <- theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  panel.grid = element_blank(),
  axis.line = element_line(color = "black", size = 0.5),
  axis.ticks = element_line(color = "black", size = 0.5),
  axis.ticks.length = unit(1.5, "mm"),
  axis.text = element_text(size = 9, color = "black"),
  axis.title = element_text(size = 11, color = "black"),
  plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
  legend.position = "none"
)

# 九个基因，按你要求的顺序
gene_list <- c("CD69",
               "CD74","HLA-DRA","HLA-DPA1","HLA-DQB1","HLA-DQA1",
               "JUN","JUNB","FOS")

# 循环生成每个基因的小提琴图
plot_list <- lapply(seq_along(gene_list), function(i) {
  g <- gene_list[i]
  
  # 根据位置决定是否显示坐标轴标题
  y_lab <- ifelse(i == 1, "Expression Level", "")
  x_lab <- "Cluster"   # 或者改成你的分组变量名，如 "Condition"
  
  p <- VlnPlot(HSC, features = g, pt.size = 0, 
               group.by = "seurat_clusters",   # 如果你的分组是二元变量，比如 "group"，改这里
               cols = group_colors) +
    # 添加小白圆点表示中位数
    stat_summary(fun = median, geom = "point", 
                 shape = 21, size = 2.5, fill = "white", color = "black", stroke = 0.4) +
    labs(title = g, x = x_lab, y = y_lab) +
    top_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  return(p)
})

# 组合成一行 9 个图
combined <- wrap_plots(plot_list, ncol = 9) +
  plot_annotation(tag_levels = 'A')  # 可选的 A-I 标签

# 保存为 PDF
ggsave("Fig8_9genes_onerow.pdf", combined, width = 18, height = 3)


# ================= 加载必要包 =================
library(Seurat)
library(GSVA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# ================= 1. 载入数据并准备表达矩阵 =================
load("D:/aJMML/scRNA2/human_HSC.RData")
HSC <- JoinLayers(HSC)   # 你已经执行过，确保可以再运行一次

# 提取 log‑normalized 数据（Seurat v5 用 layer 替代 slot）
expr <- as.matrix(GetAssayData(HSC, assay = "RNA", layer = "data"))

# ================= 2. 手动构建 IL‑10 / JAK‑STAT / 细胞因子 基因集 =================
pathways_list <- list(
  HALLMARK_IL6_JAK_STAT3_SIGNALING = c(
    "IL6ST","JAK1","JAK2","STAT3","SOCS3","IL6R","IL6","OSM","LIF",
    "CNTF","CSF2RB","CSF3R","IL10RA","IL10RB","IL4R","IL2RG",
    "STAT1","STAT5A","STAT5B","PTPN11","SOCS1","MYD88","IRAK1",
    "MAPK1","MAPK3","NFKB1","RELA","TNFRSF1A"
  ),
  HALLMARK_IL2_STAT5_SIGNALING = c(
    "IL2RA","IL2RB","IL2RG","JAK1","JAK3","STAT5A","STAT5B",
    "BCL2","BCL2L1","CCND2","CCND3","MYC","IL4R","IL7R",
    "IL15RA","IL21R","SOCS1","SOCS2","CISH","PIM1","CD69"
  ),
  HALLMARK_INFLAMMATORY_RESPONSE = c(
    "IL1B","IL1R1","TNF","TNFRSF1A","IL6","IL6R","CXCL8",
    "CCL2","CCL3","CCL5","NFKB1","RELA","STAT3","JUN","FOS",
    "MYD88","TLR2","TLR4","ICAM1","VCAM1","PTGS2","IL10"
  ),
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = c(
    "TNF","TNFRSF1A","NFKB1","RELA","NFKBIA","IKBKB",
    "JUN","JUNB","FOS","FOSL1","IL6","CXCL1","CXCL2",
    "CCL2","CCL5","ICAM1","VCAM1","BCL2A1","BCL2L1"
  ),
  KEGG_JAK_STAT_SIGNALING_PATHWAY = c(
    "IL2RA","IL2RB","IL2RG","IL3RA","IL4R","IL5RA","IL6R",
    "IL6ST","IL7R","IL9R","IL10RA","IL10RB","IL11RA","IL12RB1",
    "IL12RB2","IL13RA1","IL15RA","IL21R","IL23R","IFNAR1",
    "IFNAR2","IFNGR1","IFNGR2","CSF2RB","CSF3R","EPOR",
    "MPL","GH1","PRLR","JAK1","JAK2","JAK3","TYK2",
    "STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6",
    "SOCS1","SOCS2","SOCS3","CISH","PIM1","MYC","CCND1",
    "BCL2","BCL2L1","AKT1","PIK3CA","MAPK1","MAPK3"
  ),
  KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION = c(
    "IL1A","IL1B","IL1R1","IL1R2","IL1RAP","IL2","IL2RA",
    "IL2RB","IL2RG","IL3","IL4","IL4R","IL5","IL6","IL6R",
    "IL7","IL7R","IL9","IL10","IL10RA","IL10RB","IL11",
    "IL12A","IL12B","IL13","IL15","IL17A","IL18","IL21",
    "IL22","TNF","TNFRSF1A","LTB","IFNG","CSF2","CSF3",
    "CXCL8","CCL2","CCL5","CXCR4","CCR5"
  ),
  REACTOME_SIGNALING_BY_INTERLEUKINS = c(
    "IL1A","IL1B","IL1R1","IL2","IL2RA","IL2RB","IL2RG",
    "IL3","IL4","IL4R","IL5","IL6","IL6R","IL6ST","IL7",
    "IL7R","IL9","IL9R","IL10","IL10RA","IL10RB","IL11",
    "IL12A","IL12B","IL13","IL15","IL17A","IL18","IL21",
    "JAK1","JAK2","JAK3","TYK2","STAT1","STAT3","STAT5A",
    "STAT5B","SOCS1","SOCS3","BCL2","MYC"
  ),
  REACTOME_INTERLEUKIN_10_SIGNALING = c(
    "IL10","IL10RA","IL10RB","JAK1","TYK2","STAT3","SOCS3",
    "BCL2L1","CCND2","MYC","IL4R"
  ),
  REACTOME_JAK_STAT_SIGNALING_AFTER_INTERLEUKIN_6 = c(
    "IL6","IL6R","IL6ST","JAK1","JAK2","STAT1","STAT3",
    "SOCS3","PTPN11","MAPK1","MAPK3"
  ),
  WP_IL10_SIGNALING_PATHWAY = c(
    "IL10","IL10RA","IL10RB","JAK1","TYK2","STAT3","STAT1",
    "SOCS1","SOCS3","BCL2","BCL2L1","MYC","CCND1","CCND3"
  ),
  WP_JAK_STAT_SIGNALING_PATHWAY = c(
    "JAK1","JAK2","JAK3","TYK2","STAT1","STAT2","STAT3",
    "STAT4","STAT5A","STAT5B","STAT6","SOCS1","SOCS3","CISH",
    "PTPN11","IFNAR1","IFNAR2","IL6ST","IL2RG"
  ),
  WP_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM = c(
    "IL2","IL2RA","IL2RB","IL2RG","IL4","IL6","IL10",
    "IFNG","TNF","JAK1","JAK2","JAK3","STAT1","STAT3",
    "STAT5A","STAT5B","SOCS1","SOCS3"
  ),
  IL10_SIGNALING_CORE = c(
    "IL10RA","IL10RB","JAK1","TYK2","STAT3","SOCS3",
    "BCL2L1","IL10","STAT1","STAT5A","STAT5B"
  )
)
# 这一步能避免因基因缺失导致的错误
pathways_list_filt <- lapply(pathways_list, function(genes) {
  intersect(genes, rownames(expr))
})
# 移除基因数小于 5 的基因集
pathways_list_filt <- pathways_list_filt[lengths(pathways_list_filt) >= 5]

# ================= 3. 运行 GSVA =================
set.seed(123)
gsva_par <- gsvaParam(
  expr = expr,
  geneSets = pathways_list_filt,
  kcdf = "Gaussian",       # log‑normalized 数据用 Gaussian
  minSize = 5,
  maxSize = 500
)

# 运行 GSVA
set.seed(123)
gsva_res <- gsva(gsva_par)   # 新版直接传入 gsvaParam 对象

# 将 GSVA 评分矩阵存入 Seurat 对象
HSC[["GSVA"]] <- CreateAssayObject(data = gsva_mat)

# ================= 4. 划分 CD69high / CD69low 组 =================
cd69_expr <- FetchData(HSC, vars = "CD69")[, 1]
thr <- median(cd69_expr)   # 中位数分界，也可改为 quantile(cd69_expr, 0.75)
HSC$CD69_group <- ifelse(cd69_expr > thr, "CD69high", "CD69low")
table(HSC$CD69_group)

# ================= 5. 差异通路活性检验 =================
path_names <- rownames(gsva_mat)

diff_list <- lapply(path_names, function(p) {
  scores <- FetchData(HSC, vars = p)[, 1]
  high <- scores[HSC$CD69_group == "CD69high"]
  low  <- scores[HSC$CD69_group == "CD69low"]
  wt <- wilcox.test(high, low)
  data.frame(
    pathway    = p,
    p_val      = wt$p.value,
    mean_diff  = mean(high) - mean(low),
    mean_high  = mean(high),
    mean_low   = mean(low),
    stringsAsFactors = FALSE
  )
})
diff_df <- bind_rows(diff_list) %>%
  mutate(p_adj = p.adjust(p_val, method = "BH")) %>%
  arrange(p_adj)

# ================= 6. 火山图：高亮 IL10 / JAK / STAT 通路 =================
highlight <- diff_df %>%
  filter(grepl("IL10|IL.10|JAK|STAT|CYTOKINE|INTERLEUKIN", pathway, ignore.case = TRUE))

p_volcano <- ggplot(diff_df, aes(x = mean_diff, y = -log10(p_adj))) +
  geom_point(color = "grey80", size = 1.5) +
  geom_point(data = highlight, color = "#E41A1C", size = 2.8) +
  geom_text_repel(data = highlight, aes(label = pathway),
                  size = 3, max.overlaps = 25, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Mean activity difference (CD69high – CD69low)",
       y = expression(-log[10]~adjusted~italic(P)),
       title = "GSVA pathway activity: CD69high vs CD69low HSC") +
  theme_minimal(base_size = 12)

print(p_volcano)

# ================= 7. 核心通路箱线图（IL‑10 与 JAK‑STAT） =================
# 选择两个代表性通路
key_paths <- c("REACTOME_INTERLEUKIN_10_SIGNALING",
               "KEGG_JAK_STAT_SIGNALING_PATHWAY")

plot_data <- FetchData(HSC, vars = c(key_paths, "CD69_group")) %>%
  pivot_longer(-CD69_group, names_to = "pathway", values_to = "activity")

p_box <- ggplot(plot_data, aes(x = CD69_group, y = activity, fill = CD69_group)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  facet_wrap(~ pathway, scales = "free_y") +
  scale_fill_manual(values = c("CD69high" = "#E41A1C", "CD69low" = "#377EB8")) +
  labs(x = "", y = "GSVA enrichment score",
       title = "IL‑10 & JAK/STAT signaling in CD69high HSC") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_box)


# =============== 更严谨的效应量可视化（Bootstrap 置信区间）===============
# ============== 修正后的 Bootstrap 效应量可视化 ==============
library(boot)
library(dplyr)
library(ggplot2)

# 6 条显著通路
sig_paths <- c(
  "HALLMARK-IL2-STAT5-SIGNALING",
  "HALLMARK-TNFA-SIGNALING-VIA-NFKB",
  "KEGG-JAK-STAT-SIGNALING-PATHWAY",
  "KEGG-CYTOKINE-CYTOKINE-RECEPTOR-INTERACTION",
  "REACTOME-SIGNALING-BY-INTERLEUKINS",
  "WP-CYTOKINE-SIGNALING-IN-IMMUNE-SYSTEM"
)

# 提取评分矩阵（细胞 × 通路）和分组向量
gsva_sub <- t(gsva_scores[sig_paths, , drop = FALSE])   # 细胞为行，通路为列
group_vec <- HSC$CD69_group   # 与 gsva_sub 行顺序一致，长度相同

# 修改 Bootstrap 统计函数：data 现在是向量，用一维索引
boot_mean_diff <- function(data, indices, group_vec) {
  d <- data[indices]          # 重抽样后的通路评分
  g <- group_vec[indices]     # 对应的分组
  mean(d[g == "CD69high"]) - mean(d[g == "CD69low"])
}

# 对每条通路做 2000 次 Bootstrap
boot_results <- lapply(sig_paths, function(p) {
  vals <- gsva_sub[, p]   # 该通路在所有细胞中的评分向量
  if (any(is.na(vals) | is.infinite(vals))) return(NULL)
  
  boot_obj <- boot(data = vals, statistic = boot_mean_diff, R = 2000,
                   group_vec = group_vec)
  ci <- boot.ci(boot_obj, type = "perc")
  if (is.null(ci)) return(NULL)
  
  data.frame(
    pathway   = p,
    mean_diff = boot_obj$t0,
    lower     = ci$percent[4],
    upper     = ci$percent[5]
  )
}) %>% bind_rows()

# 绘制森林图
p_boot <- ggplot(boot_results, aes(x = mean_diff, y = reorder(pathway, mean_diff))) +
  geom_point(color = "#E41A1C", size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2,
                 color = "grey30", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Mean GSVA score difference (CD69high – CD69low)",
       y = "",
       title = "CD69high HSC pathway activation (bootstrap 95% CI)") +
  theme_bw(base_size = 11)

print(p_boot)
# 保存为 PDF（矢量，适合论文）
ggsave("D:/aJMML/scRNA2/Fig8_bootstrap_forest.pdf",
       plot = p_boot, width = 8, height = 5, device = "pdf")
