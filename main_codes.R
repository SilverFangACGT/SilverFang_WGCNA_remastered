# Setup working directory
setwd('C:/Users/Monika/Desktop/GSL')

# Read the expression data
exp <- read.table("GSL.txt", header = TRUE, sep = "\t", row.names = 1)

# Filter out low-expressed genes
gene <- subset(exp, rowSums(exp)/ncol(exp) >= 0.5)

# Transpose matrix (samples in rows, genes in columns)
gene <- t(gene)

# Perform hierarchical clustering on the samples
sampleTree <- hclust(dist(gene), method = "average")

# Plot a sample clustering tree to detect outliers
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

# Set the soft threshold range
powers <- 1:20

# Select a Network Topology
type <- "unsigned"

# Select a proper soft threshold value
sft <- pickSoftThreshold(gene, powerVector = powers, networkType = "unsigned", verbose = 5)

# Set the graphic layout to 1 row and 2 columns
par(mfrow = c(1, 2))

# 绘制尺度无关性曲线 Plotting Scale-Invariant Curves
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type = 'n', 
     xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology Model Fit,signed R^2', 
     main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, col = 'red')
abline(h = 0.8, col = 'red')

# 绘制平均连通性曲线 Plotting the Average Connectivity Curve
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n', 
     main = paste('Mean connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = 'red')

# Constructing a Gene Co-expression Network
net <- blockwiseModules(gene, power = 13, 
                        maxBlockSize = 20000, TOMType = "unsigned", 
                        minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, 
                        numericLabels = TRUE, pamRespectsDendro = FALSE, 
                        saveTOMs = FALSE, verbose = 3)

# View the number of modules
table(net$colors)

# Convert module color labels to visual colors
dynamicColors <- labels2colors(net$colors)

# View the number of dynamic color modules
table(dynamicColors)

# Plotting Gene Trees and Module Colors
plotDendroAndColors(net$dendrograms[[1]], dynamicColors[net$blockGenes[[1]]],
                    "dynamicColors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# 计算模块特征基因（ME）
MEList <- moduleEigengenes(gene, colors = dynamicColors)
MEs <- MEList$eigengenes

# 保存模块特征基因到文件
write.table(MEs, 'moduleEigengenes.txt', sep = '\t', col.names = NA, quote = FALSE)

# 计算模块特征基因之间的相关性
ME_cor <- cor(MEs)
ME_cor[1:6, 1:6]

# 对模块特征基因进行层次聚类
METree <- hclust(as.dist(1 - ME_cor), method = 'average')

# 绘制模块特征基因聚类树
plot(METree, main = 'Clustering of module eigengenes', xlab = '', sub = '')
abline(h = 0.2, col = 'blue')
abline(h = 0.25, col = 'red')

# 绘制模块特征基因网络图
plotEigengeneNetworks(MEs, '', cex.lab = 0.8, xLabelsAngle= 90,
                      marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2))

# 合并相似模块
merge_module <- mergeCloseModules(gene, dynamicColors, cutHeight = 0.2, verbose = 3)
mergedColors <- merge_module$colors

# 查看合并后的模块数量
table(mergedColors)

# 将基因和模块颜色对应关系保存到文件
gene_module <- data.frame(gene_name = colnames(gene),
                          module = mergedColors, stringsAsFactors = FALSE)
write.table(gene_module, 'gene_module.txt', sep = '\t', col.names = NA, quote = FALSE)

# 读取性状数据
trait <- read.table("Trait.txt",
                    header = TRUE,
                    row.name = 1,
                    sep = "\t")

# 获取合并后的模块特征基因
module <- merge_module$newMEs

# 计算模块与性状的相关性
moduleTraitCor <- cor(module, trait, use = 'p')

# 计算相关性显著性
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(module))

# 生成相关性矩阵的文本表示
textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep = '')
dim(textMatrix) <- dim(moduleTraitCor)

# 设置图形边距
par(mar = c(0, 0, 0, 0))

# 绘制模块与性状的热图
labeledHeatmap(Matrix = moduleTraitCor, main = paste('Module-trait'),
               xLabels = names(trait), yLabels = names(module), ySymbols = names(module),
               colorLabels = FALSE, colors = blueWhiteRed(100), cex.text = 0.6, zlim = c(-1,1),
               textMatrix = textMatrix, setStdMargins = FALSE)

# 创建用于保存Cytoscape输出的目录
dir.create('cytoscape', recursive = TRUE)

# 导出每个模块的网络到Cytoscape
for (block in 1:length(net$dendrograms)) {
  blockGenes <- net$blockGenes[[block]]
  blockColors <- dynamicColors[blockGenes]
  blockData <- gene[, blockGenes]
  blockTOM <- TOMsimilarityFromExpr(blockData, power = 13, TOMType = "unsigned")
  
  for (mod in unique(blockColors)) {
    modGenes <- which(blockColors == mod)
    modTOM <- blockTOM[modGenes, modGenes]
    dimnames(modTOM) <- list(colnames(blockData)[modGenes], colnames(blockData)[modGenes])
    outEdge <- paste('cytoscape/', mod, '_block', block, '.edge_list.txt', sep = '')
    outNode <- paste('cytoscape/', mod, '_block', block, '.node_list.txt', sep = '')
    exportNetworkToCytoscape(modTOM,
                             edgeFile = outEdge,
                             nodeFile = outNode,
                             weighted = TRUE,
                             threshold = 0.3,  
                             nodeNames = colnames(blockData)[modGenes],
                             altNodeNames = colnames(blockData)[modGenes],
                             nodeAttr = blockColors[modGenes])
  }
}
