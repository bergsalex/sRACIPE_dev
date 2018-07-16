
plot_data_knockout_single = function(output_file, output_file_knockout, plot_filename=filename,  topology_df=topology, KNOCKOUT = NA_character_, config = configuration, bin_count=40){
  message("Plotting the single knockout results")
  library(reshape2)
  library(ggplot2)
  library(Rtsne)
  library(gplots)
  library(MASS)
  library(RColorBrewer)

  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  plot_color <- rf(32)

  if(missing(bin_count)) bin_count <- 40
  if(missing(plot_filename)) output_filename <- "plots"

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


  pdf(paste("results/",plot_filename,KNOCKOUT,"_knockout.pdf",sep = ""))

  #col_count <- as.integer(65536/dim(data_simulation)[2])
  #heatmap_data <- t(data_simulation[1:col_count,])
  #dist <- as.dist((1-cor(t(heatmap_data), method = "spear"))/2)
  #cluster <- hclust(dist,method = 'ward.D2')
  #dendrogram <- as.dendrogram(cluster)
  full <- heatmap.2(t(data_simulation), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2), trace = "none")

  heatmap.2(t(data_simulation_knockout[,rev(full$rowInd)]), Rowv=NA, breaks=full$breaks, col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2), trace = "none")

  # heatmap.2(t(data_simulation_knockout), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2), trace = "none")

  dev.off()

  message("Plots in the pdf file in the results folder.")

}
