# Code to run Monocle3 through SeuratWrappers

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(scales)
library(cluster)

so = readRDS(<SeuratObject>)

# CONVERT SEURAT OBJECT TO MONOCLE3 CELL_DATA_SET OBJECT
cds = as.cell_data_set(so)

# USE MONOCLE3 TO CLUSTER AND LEARN CONNECTIVITY GRAPH
cds = cluster_cells(cds, reduction_method="UMAP") #This dimensional reduction label is case-sensitive
cds = learn_graph(cds)
plot_cells(cds)

# OPTIONS FOR PSEUDOTIME ORIGIN
# OPTION 1: USE A CELL POPULATION MATCHING A FEATURE IN THE SEURAT OBJECT, E.G. CLUSTER
originPop = colnames(so)[which(so$SCT_snn_res.1.2 == 15)]
cds = order_cells(cds, reduction_method="UMAP",root_cells = originPop)
plot_cells(cds,color_cells_by="pseudotime")

# OPTION 2: PICK A RANDOM CELL FROM THE CLUSTER TO START THE ORIGIN
originPop = colnames(so)[which(so$SCT_snn_res.1.2 == 15)]
cds = order_cells(cds, reduction_method="UMAP",root_cells = originPop[sample(length(originPop))])
plot_cells(cds,color_cells_by="pseudotime")

# OPTION 3: TRY TO USE THE NEAREST VERTEX TO THE CLUSTER AS THE ORIGIN
cell_ids = which(colData(cds)[,"SCT_snn_res.1.2"]==15)
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds),])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
cds = order_cells(cds,reduction_method="UMAP",root_pr_nodes=root_pr_nodes)

#FINISH BY ADDING PSEUDOTIME TO SEURAT OBJECT AS METADATA
so = AddMetaData(so,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name="pseudotime"
)
FeaturePlot(so,features=mesoderm) & scale_color_viridis_c()
