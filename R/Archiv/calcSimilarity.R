
calcSimilarity = function(data_simulation, data_simulation_knockout, n_clusters){
  message("Finding the similarity index")
  library(reshape2)
  library(ggplot2)
  library(Rtsne)
  library(gplots)
  library(MASS)
  library(RColorBrewer)
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  plot_color <- rf(32)

  rm(list=ls())
  topology_file <- "/Users/koharv/Documents/Work/ScenicMay/IPSC/output/IPSC_network_w001_corr002_Target.txt"
  topology <- sRACIPE_load_topology(topology_file)
  KNOCKOUT = "ATF3"
  output_file <- "/Users/koharv/Documents/GitHub/sRACIPE_dev/results/sRACIPE_RK_IPSC_network_w001_corr002_Target_g31_output.txt"
  output_file_knockout <- paste("/Users/koharv/Documents/GitHub/sRACIPE_dev/results/sRACIPE_RK_IPSC_network_w001_corr002_Target_",KNOCKOUT,"_KO_g31_output.txt",sep = "")
  working_directory <- getwd()
  name_genes <- read.table(paste(working_directory,"/results/gene_interaction_topology_",topology$filename,".txt",sep=""), header = T, stringsAsFactors = F)
  data_simulation <- read.table(output_file, header = F)
  data_simulation <- load_data(data_simulation, topology)

  name_genes <- t(as.matrix(name_genes))

  knockout_number <- as.integer(which(name_genes==KNOCKOUT))
  data_simulation <- data_simulation[,-knockout_number]

  data_simulation_knockout <- read.table(output_file_knockout, header = F)
  data_simulation_knockout <- load_data(data_simulation_knockout, topology)
  data_simulation_knockout <- data_simulation_knockout[,-knockout_number]


  n_models <- 5000
  n_modelsKO <- 5000
  data_simulation1 <- t(data_simulation[1:n_models,])
  data_simulation_knockout1 <- t(data_simulation_knockout[1:n_modelsKO,])

  distance <- as.dist((1-cor((data_simulation1), method = "spearman"))/2)
  clusters <- hclust(distance, method = 'ward.D2')
  #plot(clusters)
  n_clusters <- 2
  clusterCut <- cutree(clusters, n_clusters)

  clusterFreq <- table(clusterCut)/n_models
  corSpear <- t(cor(data_simulation1, data_simulation_knockout1, method = "p"))
  if(sum(is.na(corSpear)) > 0) message("Error in correlation. Please verify the data")
  clusterCor <- matrix(0, nrow=n_modelsKO, ncol = n_clusters)


  for(i in 1:n_modelsKO){
  for(j in 1:n_clusters)
  {
    clusterCor[i,j] <- mean(corSpear[i, which(clusterCut==j)])
  }
  }
  clusterKO <- (apply(clusterCor,1,which.max))
  clusterKOFreq <- table(clusterKO)/n_modelsKO

  clusterKOcor <- c(rep(0, length(clusterKO)))
  for(i in 1: nrow(clusterCor))
  {
    clusterKOcor[i] <- clusterCor[i,clusterKO[i]]
  }

  clusterKOcorMean <- c(rep(0, n_clusters))
  for(i in 1: n_clusters)
  {
    clusterKOcorMean[i] <- mean(clusterKOcor[which(clusterKO == i)])
  }

similarity <- clusterKOFreq*clusterFreq*clusterKOcorMean
sum(similarity)*sqrt(n_clusters)
#return(similarity)
}
