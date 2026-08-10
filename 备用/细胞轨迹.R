y####1. 数据预处理####

#1.1 安装monocle2
install.packages("BiocManager")
BiocManager::install("monocle")
install.packages("RColorBrewer") #画图配色板

#1.2 加载安装包
library(monocle)
library(RColorBrewer)
setwd("工作路径")
#1.3 读取数据GSE90047
# TPM矩阵（行为基因名，列为细胞名）
expr_matrix <- read.table("0.data/scRNA-seq_TPM_GSE90047.xls",header=T,row.names=1, sep = "\t")
# 细胞注释矩阵（行为细胞名）
sample_sheet <- read.table("0.data/XCR_paper_cell_anno_sub.txt",header=T,row.names=1, sep ="\t") 
# 基因注释矩阵（行为基因名）
gene_annotation <- read.table("0.data/Mus_ref.txt",header=T,row.names=1, sep = "\t")
# 从注释信息中筛选出有表达的基因
pos <- which(rownames(gene_annotation) %in% rownames(expr_matrix))
gene_annotation <- gene_annotation[pos,]
# 注意：expr_matrix中列名的顺序和 sample_sheet中行名的顺序必须一致
expr_matrix <- expr_matrix[rownames(gene_annotation),rownames(sample_sheet)]

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

####2. 创建对象####
#2.1 创建monocle2的对象
# 创建monocle对象
cd <- newCellDataSet(as.matrix(expr_matrix), 
                     phenoData = pd, 
                     featureData = fd, 
                     lowerDetectionLimit = 0.1, 
                     expressionFamily = tobit(Lower = 0.1))
# expressionFamily: 数据为TPM/FPKM时设置为tobit(Lower = 0.1)，数据为count时设置为negbinomial.size())

# 将FPKM/TPM数据转换为UMI数据（read count）
rpc_matrix <- relative2abs(cd)

# 重新创建monocle对象
cd <- newCellDataSet(as(as.matrix(rpc_matrix),"sparseMatrix"), 
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())
#2.2 计算数据的size factors和dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)
#2.3 数据过滤
# 过滤低于1%细胞中检出的基因，最低表达阈值为0.5
cd <- detectGenes(cd, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))
# 查看筛选后的基因个数（用于后面基因的筛选）
length(expressed_genes)

####3. 构建轨迹####
#3.1 关键基因的筛选
#3.1.1 挑选细胞间高度离散的基因
disp_table <- dispersionTable(cd)
# 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cd <- setOrderingFilter(cd, ordering_genes)
plot_ordering_genes(cd)
#3.1.2 挑选在指定类型间差异表达的基因
# 根据数据的某种分类（如本例中根据Stage的类型）进行差异分析
diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = "~Stage")
# 筛选差异基因（q < 1e-5并且属于之前计算的expressed_genes列表中）
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-5))
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)
plot_ordering_genes(cd)
#3.1.3 根据已知标记基因对细胞进行排序
#3.2 细胞降维
# 选取的基因数目为每个细胞的维度，基于默认的'DDRTree'方法进行数据降维
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
#3.3 细胞排序
# 对细胞进行排序，由于排序无法区分起点和终点，若分析所得时序与实际相反，根据“reverse”参数进行调整，默认reverse=F
cd <- orderCells(cd, reverse = F) 
#3.4 轨迹可视化
# 排序好的细胞可以进行可视化，可标注细胞的各注释信息"Stage", "Putative_cell_type", "Sort_by", "State", "Pseudotime"等)
# 查看细胞注释信息
head(cd@phenoData@data)
# 设置颜色
color1 <- c(brewer.pal(6, "Set1"))
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))
# 按照时间点进行映射
plot_cell_trajectory(cd, show_cell_names = F, color_by = "Stage") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"Stage"]))))
# 按照细胞类型进行映射（手动设置颜色）
plot_cell_trajectory(cd, show_cell_names = F, color_by = "Putative_cell_type") + scale_color_manual(values = c("#E41A1C", "#999999"))
# 按照State进行映射
plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "State") +  scale_color_manual(values = color1))
# 按照计算的Pseudotime进行映射
plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "Pseudotime") + scale_color_viridis_c())

