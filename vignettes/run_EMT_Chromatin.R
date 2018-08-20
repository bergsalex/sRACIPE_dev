rm(list = ls())

#plot_data_stochastic(output_file, plot_filename=topology$filename,  topology_df=topology, config = configuration, bin_count=40)


#output_file <- sRACIPE_stochastic(topology_file = topology_file, NUM_MODELS = 10000, ANNEAL = F, NOISE_LEVELS = 30, MAX_NOISE = 50,NOISE_SCALING_FACTOR = 0.5)

#output_file <- sRACIPE_stochastic(topology_file = topology_file, NUM_MODELS = 10000, ANNEAL = F, NOISE_LEVELS = 30, NOISE_SCALING_FACTOR = 0.5)

##########################################################################################################################################
# Analyze the Epcam- network

#########################################################################################################################################
library(reshape2)
library(ggplot2)
library(Rtsne)
library(gplots)
library(MASS)
library(RColorBrewer)
library(htmlwidgets)
library(d3heatmap)
library(sRACIPEv03)

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot_color <- rf(32)
setwd("/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/")


topology_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/inputs/EMT_chromatin.tpo"
topology <- sRACIPE_load_topology(topology_file = topology_file)

output_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/results/sRACIPE_EM_EMT_chromatin_g15_Annealing_0_output.txt"
data_simulation_EpcamN <- (as.data.frame(read.table(output_file, header = F)))
data_simulation_test <- data_simulation_EpcamN[,436:450]
data_simulation_test <- log2(data_simulation_test)
col.means.EpcamN <- colMeans(data_simulation_test)
col.sds.EpcamN <- apply(data_simulation_test, 2, sd)
name_genes.EpcamN <- read.table(paste(getwd(),"/results/gene_interaction_topology_",topology$filename,".txt",sep=""), header = T, stringsAsFactors = F)
name_genes.EpcamN <- t(as.matrix(name_genes.EpcamN))

data_simulation_test <- data_simulation_EpcamN[1:5000,436:450]
name_models.EpcamN <- seq(1:nrow(data_simulation_test))
row.names(data_simulation_test) <- name_models.EpcamN
colnames(data_simulation_test) <- name_genes.EpcamN

data_simulation_test <- log2(data_simulation_test)
data_simulation_test <- data_simulation_test[is.finite(rowMeans(data_simulation_test)), ]

data_simulation_test <- sweep(data_simulation_test,2,col.means.EpcamN,FUN = "-")
data_simulation_test <- sweep(data_simulation_test,2,col.sds.EpcamN,FUN = "/")

heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'complete'))

#heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2))

#####################################
## EpcamN with annealing
output_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/results/sRACIPE_EM_EMT_chromatin_g15_Annealing_1_output.txt"
data_simulation_EpcamN_A <- (as.data.frame(read.table(output_file, header = F)))

data_simulation_test <- data_simulation_EpcamN_A[1:5000,436:450]
name_models.EpcamN <- seq(1:nrow(data_simulation_test))
row.names(data_simulation_test) <- name_models.EpcamN
colnames(data_simulation_test) <- name_genes.EpcamN

data_simulation_test <- log2(data_simulation_test)
data_simulation_test <- data_simulation_test[is.finite(rowMeans(data_simulation_test)), ]

data_simulation_test <- sweep(data_simulation_test,2,col.means.EpcamN,FUN = "-")
data_simulation_test <- sweep(data_simulation_test,2,col.sds.EpcamN,FUN = "/")

heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'complete'))
##########################################################################################################################################
# Analyze the Epcam+ network

#########################################################################################################################################

topology_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/inputs/EMT_chromatinEpcamP.tpo"
topology_EpcamP <- sRACIPE_load_topology(topology_file = topology_file)

output_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/results/sRACIPE_EM_EMT_chromatinEpcamP_g11_Annealing_0_output.txt"
data_simulation_EpcamP <- (as.data.frame(read.table(output_file, header = F)))
data_simulation_test <- data_simulation_EpcamP[,320:330]
data_simulation_test <- log2(data_simulation_test)
col.means.EpcamP <- colMeans(data_simulation_test)
col.sds.EpcamP <- apply(data_simulation_test, 2, sd)
name_genes.EpcamP <- read.table(paste(getwd(),"/results/gene_interaction_topology_",topology_EpcamP$filename,".txt",sep=""), header = T, stringsAsFactors = F)
name_genes.EpcamP <- t(as.matrix(name_genes.EpcamP))

data_simulation_test <- data_simulation_EpcamP[1:5000,320:330]
name_models.EpcamP <- seq(1:nrow(data_simulation_test))
row.names(data_simulation_test) <- name_models.EpcamP
colnames(data_simulation_test) <- name_genes.EpcamP

data_simulation_test <- log2(data_simulation_test)
data_simulation_test <- data_simulation_test[is.finite(rowMeans(data_simulation_test)), ]

data_simulation_test <- sweep(data_simulation_test,2,col.means.EpcamP,FUN = "-")
data_simulation_test <- sweep(data_simulation_test,2,col.sds.EpcamP,FUN = "/")

heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'complete'))

#heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2))

#####################################
## EpcamP with annealing
output_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/results/sRACIPE_EM_EMT_chromatinEpcamP_g11_Annealing_1_output.txt"
data_simulation_EpcamP_A <- (as.data.frame(read.table(output_file, header = F)))

data_simulation_test <- data_simulation_EpcamP_A[1:5000,320:330]
name_models.EpcamP <- seq(1:nrow(data_simulation_test))
row.names(data_simulation_test) <- name_models.EpcamP
colnames(data_simulation_test) <- name_genes.EpcamP

data_simulation_test <- log2(data_simulation_test)
data_simulation_test <- data_simulation_test[is.finite(rowMeans(data_simulation_test)), ]

data_simulation_test <- sweep(data_simulation_test,2,col.means.EpcamP,FUN = "-")
data_simulation_test <- sweep(data_simulation_test,2,col.sds.EpcamP,FUN = "/")

heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'complete'))
##########################################################################################################################################
# Analyze the Epcam All network

#########################################################################################################################################

topology_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/inputs/EMT_chromatin_all.tpo"
topology_all <- sRACIPE_load_topology(topology_file = topology_file)

output_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/results/sRACIPE_EM_EMT_chromatin_all_g16_Annealing_0_pmFile_0_output.txt"
data_simulation_All <- (as.data.frame(read.table(output_file, header = F)))
data_simulation_test <- data_simulation_All[,145:160]
data_simulation_test <- log2(data_simulation_test)
col.means.All <- colMeans(data_simulation_test)
col.sds.All <- apply(data_simulation_test, 2, sd)
name_genes.All <- read.table(paste(getwd(),"/results/gene_interaction_topology_",topology_all$filename,".txt",sep=""), header = T, stringsAsFactors = F)
name_genes.All <- t(as.matrix(name_genes.All))

data_simulation_test <- data_simulation_All[1:2000,97:112]
name_models.All <- seq(1:nrow(data_simulation_test))
row.names(data_simulation_test) <- name_models.All
colnames(data_simulation_test) <- name_genes.All

data_simulation_test <- log2(data_simulation_test)
data_simulation_test <- data_simulation_test[is.finite(rowMeans(data_simulation_test)), ]

data_simulation_test <- sweep(data_simulation_test,2,col.means.All,FUN = "-")
data_simulation_test <- sweep(data_simulation_test,2,col.sds.All,FUN = "/")

heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'))

heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2))


#pca1 = prcomp((data_simulation_test), scale. =T)
#pca_data <- data.frame(x=pca1$x[,1],y=pca1$x[,2])

pca_data <- scale(data_simulation_test, pca1$center, pca1$scale) %*% pca1$rotation
pca_data <- data.frame(x=pca_data[,1],y=pca_data[,2])

pca_data$z <- data_simulation_test$Epcam
ggplot(pca_data, aes(x=x, y=y ) ) +
  #scale_fill_gradient(low="blue", high="red") +
  geom_point(aes(color=z), size=1) +
  scale_color_gradient(low = "darkblue", high = "orange",  limits=c(-1.5,1.5)) +
#scale_fill_distiller(palette= "Spectral", direction = -1) +

  #stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #geom_contour((aes(x=x, y=y,z = ..density..))) +
  #scale_fill_distiller(palette= "Spectral", direction = -1)  +
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  #xlim(-15,20) +
  #ylim(-15,20) +
  labs(x = 'Normalized Y expression('~log[2]~')', y = 'Normalized X expression('~log[2]~')') +
  theme(
    legend.position='none',
    text = element_text(size=20),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


#####################################
## Epcam All with annealing
output_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/results/sRACIPE_EM_EMT_chromatin_all_g16_Annealing_1_pmFile_0_output.txt"
data_simulation_All_A <- (as.data.frame(read.table(output_file, header = F)))

data_simulation_test <- data_simulation_All_A[1:5000,465:480]
name_models.All <- seq(1:nrow(data_simulation_test))
row.names(data_simulation_test) <- name_models.All
colnames(data_simulation_test) <- name_genes.All

data_simulation_test <- log2(data_simulation_test)
data_simulation_test <- data_simulation_test[is.finite(rowMeans(data_simulation_test)), ]

data_simulation_test <- sweep(data_simulation_test,2,col.means.All,FUN = "-")
data_simulation_test <- sweep(data_simulation_test,2,col.sds.All,FUN = "/")

#heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'))
heatmap(t(data_simulation_test), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2))

###############################################################################################################
###############################################################################################################

output_file <- "/Users/koharv/Documents/Work/sRACIPE_project/EMT_nature/results/sRACIPE_EM_EMT_chromatin_all_g16_Annealing_0_pmFile_0_output.txt"


test <- data_simulation_test


bin_count = 40
plot_data <- data.frame(x=test[,1],y=test[,2])
colnames(plot_data) <- c("x", "y")

ggplot(plot_data, aes(x=y, y=x) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #geom_contour((aes(x=x, y=y,z = ..density..))) +
  scale_fill_distiller(palette= "Spectral", direction = -1)  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #  xlim(-1.5,1.5) +
  #  ylim(-1.5,1.5) +
  labs(x = 'Normalized Y expression('~log[2]~')', y = 'Normalized X expression('~log[2]~')') +
  theme(
    legend.position='none',
    text = element_text(size=20)
  )


 bin_count <- 40
output_filename <- "plots"

working_directory <- getwd()


data_simulation_all <- (as.data.frame(read.table(output_file, header = F)))
test <- as.matrix(data_simulation_all, nrow = 10000, ncol = 450)
dim(test) <- c(30000, 15)
d <- t(d)

col_start <- topology$number_gene*(configuration$NOISE_LEVELS-1)+1
col_end <- topology$number_gene*(configuration$NOISE_LEVELS)
data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
data_simulation <- log2(data_simulation)
data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]

means <- colMeans(data_simulation)
sds <- apply(data_simulation, 2, sd)
data_simulation <- sweep(data_simulation,2,means,FUN = "-")
data_simulation <- sweep(data_simulation,2,sds,FUN = "/")

name_models <- seq(1:nrow(data_simulation))
name_genes <- t(as.matrix(name_genes))
row.names(data_simulation) <- name_models
colnames(data_simulation) <- name_genes

i=1
col_start <- (topology$number_gene)*(configuration$NOISE_LEVELS-i)+1
col_end <- (topology$number_gene)*(configuration$NOISE_LEVELS-i+1)
data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
data_simulation <- log2(data_simulation)
data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]

data_simulation <- sweep(data_simulation,2,means,FUN = "-")
data_simulation <- sweep(data_simulation,2,sds,FUN = "/")

name_models <- seq(1:nrow(data_simulation))
row.names(data_simulation) <- name_models
colnames(data_simulation) <- name_genes

#heatmap(t(data_simulation), col=plot_color, hclustfun = function(x) hclust(x,method = 'complete'))

heatmap(t(data_simulation), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2))


pca = prcomp(data_simulation, scale. = FALSE, center = FALSE)
#pca_data <- data.frame(x=pca$x[,1],y=pca$x[,2])

pdf(paste("results/",plot_filename,"_stochastic_pca.pdf",sep = ""))

for(i in 1:configuration$NOISE_LEVELS)
{
  col_start <- (topology$number_gene)*(configuration$NOISE_LEVELS-i)+1
  col_end <- (topology$number_gene)*(configuration$NOISE_LEVELS-i+1)
  data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
  data_simulation <- log2(data_simulation)
  data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]

  data_simulation <- sweep(data_simulation,2,means,FUN = "-")
  data_simulation <- sweep(data_simulation,2,sds,FUN = "/")

  name_models <- seq(1:nrow(data_simulation))
  row.names(data_simulation) <- name_models
  colnames(data_simulation) <- name_genes

  pca_data <- scale(data_simulation, pca$center, pca$scale) %*% pca$rotation
  pca_data <- data.frame(x=pca_data[,1],y=pca_data[,2])
  density_plot(pca_data,bin_count,plot_color)
}


dev.off()





plot_data_deterministic(output_file, plot_filename = topology$filename)
