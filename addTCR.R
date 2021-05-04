library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
rds <- as.character(args[1])
vdj <- as.character(args[2])

pbmc = readRDS(paste0("rds/",rds))
meta = read.delim(vdj,sep = ',',stringsAsFactors = F)
meta$barcode = gsub("-1","",meta$barcode)
meta$cdr3_length = nchar(meta$cdr3)
meta$cdr3_nt_length = nchar(meta$cdr3_nt)

df.agg <- aggregate(reads ~ barcode+chain, meta, max)
meta <- merge(df.agg, meta)

tra=meta[meta$chain == "TRA",]
trb=meta[meta$chain == "TRB",]

newMeta = merge(tra,trb,by="barcode",all=T)
rownames(newMeta) = make.names(newMeta$barcode,unique = T)

newMeta$matching.x = "No match"
for (tcrX in newMeta$cdr3.x) {
  counts = head(sort(table(newMeta$cdr3.y[which(newMeta$cdr3.x == tcrX)]),decreasing = T),1)
  occur=paste(names(counts),counts)
  if (length(occur) > 0) 
  newMeta$matching.x[which(newMeta$cdr3.x == tcrX)] = occur
}

newMeta$matching.y = "No match"
for (tcrY in newMeta$cdr3.y) {
  counts = head(sort(table(newMeta$cdr3.x[which(newMeta$cdr3.y == tcrY)]),decreasing = T),1)
  occur=paste(names(counts),counts)
  if (length(occur) > 0) 
    newMeta$matching.y[which(newMeta$cdr3.y == tcrY)] = occur
}

rownames(newMeta) = paste0(rownames(newMeta),"-1")

pbmc = AddMetaData(pbmc, newMeta)

saveRDS(pbmc, paste0("tcr/",rds))
