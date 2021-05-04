library(Seurat)
library(slingshot)

runSlingShot = function(obj,seperateBy,output) {
  require(slingshot)
  require(RColorBrewer)
  require(SingleCellExperiment)
  require(gam)
  require(clusterExperiment)
  require(Seurat)
  sim = as.SingleCellExperiment(obj)
  sce <- slingshot(sim, clusterLabels = seperateBy, reducedDim = 'PCA')
  pdf(paste0(output,"_slingshot.pdf"))
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
  plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, col='black')
  plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
  t <- sce$slingPseudotime_1
  # for time, only look at the 100 most variable genes
  Y <- log1p(assays(sim)$logcounts)
  var1K <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
  Y <- Y[var1K,]
  # fit a GAM with a loess term for pseudotime
  gam.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    tmp <- gam(z ~ lo(t), data=d)
    p <- summary(tmp)[4][[1]][1,5]
    p
  })
  topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
  heatdata <- assays(sim)$logcounts[rownames(assays(sim)$logcounts) %in% topgenes,
                                    order(t, na.last = NA)]
  heatclus <- sce[[seperateBy]][order(t, na.last = NA)]
  ce <- ClusterExperiment(as.matrix(heatdata), heatclus, transformation = log1p)
  plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
              visualizeData = 'transformed')
  rd=obj@reductions$umap@cell.embeddings[,1:2]
  cl = sim[[seperateBy]]
  lin1 <- getLineages(rd, cl)
  plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
  lines(lin1, lwd = 3, col = 'black')
  crv1 <- getCurves(lin1)
  plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
  lines(crv1, lwd = 3, col = 'black')
  dev.off()
}
load(".RData")
runSlingShot(ccbr970_newDat.merged,"seurat_clusters","slingshot_ccbr970")

