library(clusterSim)

embeddingMat = as.data.frame(Embeddings(so)[,1:30])
distance = dist(embeddingMat)
clusters = as.numeric(as.character(so$SCT_snn_res.0.8))

#Calinski-Harabasz index: Higher scores preferred
ch_score=clusterSim::index.G1(x=embeddingMat,cl=clusters)

#Silhouette score: Higher scores preferred
sil_score = clusterSim::index.S(d=distance,cl=clusters)

#Davies-Bouldin index: Lower scores preferred
db_score = clusterSim::index.DB(x=embeddingMat,cl=clusters)$DB
