import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import loompy
import os 
from matplotlib.backends.backend_pdf import PdfPages
#from PyPDF2 import PdfFileMerger
import shutil
#import scvelo as scv


### EXPORT DATA FROM R Seurat objects (assuming that all normalization, dimensionality reduction, clustering has been done already)
#library(Seurat)
#library(Matrix)

## EXPORTING GENE DATA
## As raw counts:
#write.table(t(as.matrix(so@assays$RNA@counts)),file = "allExpression.txt", sep="\t",quote=F)#Note: if SCANPY can import sparse matrices, leave out the as.matrix call. Rownames should be cell IDs
## As log normalized counts 
#write.table(t(as.matrix(so@assays$RNA@data)),file = "allExpression.txt", sep="\t",quote=F)
## As SCT normalized counts
#write.table(t(as.matrix(so@assays$SCT@data)),file = "allExpression.txt", sep="\t",quote=F)

## EXPORTING UMAP COORDINATES
#write.table(so@reductions$umap@cell.embeddings, file="umap_coords.txt",sep="\t",quote=F)

## EXPORTING INITIAL CLUSTERS
## Note: You can apply this to any metadata in the Seurat object
#write.table(so$SCT_snn_res.0.6,file = "allClust.txt",sep = "\t", quote=F,col.names="seurat_clust")

## FINAL NOTE: You may have to edit the files directly to add a leading tab. I normally use sed in bash: `sed -i s"/^ColumnHeader/\tColumnHeader/"`

#SET UP THE ANNDATA OBJECT
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()
z = sc.read_csv("allExpression.txt",delimiter="\t")
sc.AnnData(z)
sc.pp.scale(z)
#sc.pp.filter_genes_dispersion(z, flavor='cell_ranger', n_top_genes=1000, log=False)
sc.tl.pca(z, svd_solver='auto')
sc.pp.neighbors(z, n_neighbors=30, n_pcs=15)
#sc.tl.louvain(z)
#sc.tl.draw_graph(z)

#IMPORT CLUSTERS
clust = pd.read_csv("allClust.txt",delimiter="\t",index_col= 0)
z.obs=pd.merge(clust,z.obs,left_index=True, right_index=True)
z.obs['seurat_clust'] = z.obs['seurat_clust'].apply(str)
#sc.pl.draw_graph(z, color='seurat_clust', legend_loc='on data')

#RUN INITAL PAGA, based on connectivity
#If you draw the graphs without the "save" flag, it prints it to visualization (e.g. X11)
sc.tl.paga(z, groups='seurat_clust')
sc.pl.paga(z,color='seurat_clust',save="outputFileHeader")

# This is a snippet for creating a new AnnData observation based on renaming, similar to the metadata columns in Seurat. Each array name in z.obs is essentially the metadata header
z.obs['seurat_anno'] = z.obs['seurat_clust']
z.obs['seurat_anno'].cat.categories = ['0_Cd69-', '1_PreTreg', '2_Cd69+', '3_MatureTreg', '4_AltPreTreg', '5_Il2', '6_Il2', '7_AltPreTreg', '8_MatureTreg', '9_Mix_Tgd', '10_MatureTreg', '11_Cd69+/PreTreg', '12_Cd69-/Cd69+','13_Cd69+', '14_Cd69-', '15_AltPreTreg','16_AltPreTreg', '17_Il2', '18_AltPreTreg', '19_AltPreTreg']
sc.tl.paga(z,groups='seurat_anno')
sc.pl.paga(z,color='seurat_anno',threshold=0.5, fontsize=8, edge_width_scale=0.6)

#This section imported additional metadata files, and created new AnnData observations called 'subclusters'  and 'CellType'
subclusters=pd.read_csv("subclusters.txt",delimiter="\t",index_col=0)
z.obs=pd.merge(subclusters,z.obs,left_index=True, right_index=True)
z.obs['subclusters'] = z.obs['subclusters'].apply(str)

sc.tl.paga(z,groups='subclusters')
sc.pl.paga(z,color='subclusters',edge_width_scale=0.4,threshold=0.8)
sc.pl.paga(z,color='seurat_clust',edge_width_scale=0.4,threshold=0.8)

cellType=pd.read_csv("cellHashingCellType.txt",delimiter="\t",index_col=0)
z.obs=pd.merge(cellType,z.obs,left_index=True,right_index=True)
z.obs['CellType'] = z.obs['CellType'].apply(str)
sc.tl.paga(z,groups="CellType")

#INITIAL SETUP FOR PSEUDOTIME ANALYSIS, USING UMAP COORDINATES
z.uns['draw_graph']={'params': {'layout': 'fr', 'random_state': 0}}
#IMPORT COORDINATES
umap = pd.read_csv("umap_coords.txt",delimiter="\t",index_col= 0)
umap = umap[['UMAP_1','UMAP_2']]
umap = umap.rename_axis('ID').values
z.obsm['X_draw_graph_fr'] = umap
sc.pl.draw_graph(z,color='seurat_clust',legend_loc="on data")
#THIS SETS THE STARTING POINT FOR PSEUDOTIME ANALYSIS. IT SELECTS THE CELL AT THE INDICATED INDEX WITHIN THE OBSERVATION (e.g. cell index 0 for Seurat cluster 2)
z.uns['iroot']=np.flatnonzero(z.obs['seurat_clust']=='2')[0]
#THIS IS THE PSEUDOTIME CALL
sc.tl.dpt(z)
sc.pl.draw_graph(z, color='dpt_pseudotime')

# TO GENERATE SUBSETS, IF NECESSARY. NOTE RE-RUNNING ALL PREPROCESSING AFTER REMOVING SUBSETS
subset = z
subset = subset[subset.obs['seurat_clust']!='14']
subset = subset[subset.obs['seurat_clust']!='15']
subset = subset[subset.obs['seurat_clust']!='16']
subset = subset[subset.obs['seurat_clust']!='18']
subset = subset[subset.obs['seurat_clust']!='19']
# 
sc.pp.filter_genes(subset,min_cells=1)
sc.pp.scale(subset)
sc.tl.pca(subset,svd_solver='auto')
sc.pp.neighbors(subset, n_neighbors=30, n_pcs=15)
sc.tl.paga(subset,groups='seurat_anno')
sc.pl.paga(subset,color='seurat_anno',threshold=0.3,save='outliersRemoved')
#

##THESE SNIPPETS WERE USED TO TRY DIFFERENT STARTING POINTS in PSEUDOTIME ANALYSIS
# subset.uns['iroot']=np.flatnonzero(subset.obs['seurat_clust']=='0')[2000]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset,color='dpt_pseudotime',save="_pseudotime_00_with_9")
# 
# subset.uns['iroot']=np.flatnonzero(subset.obs['seurat_clust']=='1')[1900]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset,color='dpt_pseudotime',save="_pseudotime_01_with_9")
# 
# subset.uns['iroot']=np.flatnonzero(subset.obs['seurat_clust']=='2')[1500]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset,color='dpt_pseudotime',save="_pseudotime_02_with_9")
# 
# subset.uns['iroot']=np.flatnonzero(subset.obs['seurat_clust']=='3')[1500]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset,color='dpt_pseudotime',save="_pseudotime_03_with_9")
# 
# subset.uns['iroot']=np.flatnonzero(subset.obs['seurat_clust']=='4')[1000]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset,color='dpt_pseudotime',save="_pseudotime_04_with_9")
# 
# subset.uns['iroot']=np.flatnonzero(subset.obs['seurat_clust']=='5')[1000]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset,color='dpt_pseudotime',save="_pseudotime_05_with_9")
# 
# subset.uns['iroot']=np.flatnonzero(subset.obs['seurat_clust']=='6')[800]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset,color='dpt_pseudotime',save="_pseudotime_06_with_9")
# 
# subset.uns['iroot'] = np.flatnonzero(z.obs['seurat_clust']=='7')[750]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset, color='dpt_pseudotime',save="_pseudotime_07_with_9")
# 
# subset.uns['iroot'] = np.flatnonzero(z.obs['seurat_clust']=='8')[500]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset, color='dpt_pseudotime',save="_pseudotime_08_with_9")
# 
# subset.uns['iroot'] = np.flatnonzero(z.obs['seurat_clust']=='9')[500]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset, color='dpt_pseudotime',save="_pseudotime_09_with_9")
# 
# subset.uns['iroot'] = np.flatnonzero(z.obs['seurat_clust']=='10')[390]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset, color='dpt_pseudotime',save="_pseudotime_10_with_9")
# 
# subset.uns['iroot'] = np.flatnonzero(z.obs['seurat_clust']=='11')[300]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset, color='dpt_pseudotime',save="_pseudotime_11_with_9")
# 
# subset.uns['iroot'] = np.flatnonzero(z.obs['seurat_clust']=='12')[250]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset, color='dpt_pseudotime',save="_pseudotime_12_with_9")
# 
# subset.uns['iroot'] = np.flatnonzero(z.obs['seurat_clust']=='13')[190]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset, color='dpt_pseudotime',save="_pseudotime_13_with_9")
# 
# 
# subset.uns['iroot']=np.flatnonzero(subset.obs['seurat_clust']=='17')[45]
# sc.tl.dpt(subset)
# sc.pl.draw_graph(subset,color='dpt_pseudotime',save="_pseudotime_17_with_9")
# 
# 