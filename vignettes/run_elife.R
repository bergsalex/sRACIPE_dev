rm(list = ls())

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


working_directory <- getwd()
config_file = "/Users/koharv/Documents/GitHub/sRACIPE_dev_JAX/inputs/sRACIPE.cfg"
configuration <- sRACIPE_load_configuration(config_file = config_file)

topology_file <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/inputs/g2_ABact.tpo"
topology <- sRACIPE_load_topology(topology_file = topology_file)

output_file_base <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_RK_g2_ABact_g2_output.txt"
name_genes <- read.table(paste("/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/gene_interaction_topology_",topology$filename,".txt",sep=""), header = T, stringsAsFactors = F)
data_simulation <- read.table(output_file_base, header = F)
#data_simulation <- data_simulation[1:50000,]
data_simulation <- log2(data_simulation)
name_genes <- t(as.matrix(name_genes))
col.means <- colMeans(data_simulation)
col.sds <- apply(data_simulation, 2, sd)
data_simulation <- sweep(data_simulation,2,col.means,FUN = "-")
data_simulation <- sweep(data_simulation,2,col.sds,FUN = "/")


bin_count=40

plot_data <- data.frame(x=data_simulation[,1],y=data_simulation[,2])
colnames(plot_data) <- c("x", "y")
h1 <- hist(plot_data$x, breaks=bin_count, plot=F)
h2 <- hist(plot_data$y, breaks=bin_count, plot=F)
top <- max(h1$counts, h2$counts)
k <- kde2d(plot_data$x, plot_data$y, n=bin_count)

# margins
oldpar <- par()
par(mar=c(3,3,1,1))
layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
image(k, col=plot_color) #plot the image
par(mar=c(0,2,1,0))
barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
par(mar=c(2,0,0.5,1))
barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)




low_high_boundary = 0.250

parameter.file <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_RK_g2_ABact_g2_parameters.txt"
parameters <- read.table(parameter.file, header = F)

parameter.selected <- parameters[,c(13,14,15,16)]
row.names(parameter.selected) <- seq(1:1000000)

plot_data <- data.frame(x=parameter.selected[,4] ,y=data_simulation[,1])
density_plot(plot_data,bin_count,plot_color)


parameter.selected.LL <- parameter.selected[which(data_simulation[,1] < low_high_boundary),]
data_simulation.LL <- data_simulation[which(data_simulation[,1] < low_high_boundary),]

parameter.selected.LL <- parameter.selected.LL[which(data_simulation.LL[,2] < low_high_boundary),]
data_simulation.LL <- data_simulation.LL[which(data_simulation.LL[,2] < low_high_boundary),]

models.LL <- rownames(parameter.selected.LL)

plot.data <- as.data.frame(parameter.selected.LL)
names(plot.data) <- c("L11", "L21", "L12", "L22")
plot.data.melted <- melt(plot.data)

ggplot(plot.data.melted, aes(x=value, color=variable)) +
  #geom_histogram(aes(y=..density..), fill="white", bins = 25, alpha=0.5, position="identity") +
  geom_density(alpha=.2)

parameter.selected.HH <- parameter.selected[which(data_simulation[,1] > low_high_boundary),]
data_simulation.HH <- data_simulation[which(data_simulation[,1] > low_high_boundary),]

parameter.selected.HH <- parameter.selected.HH[which(data_simulation.HH[,2] > low_high_boundary),]
data_simulation.HH <- data_simulation.HH[which(data_simulation.HH[,2] > low_high_boundary),]

models.HH <- rownames(parameter.selected.HH)

plot.data <- as.data.frame(parameter.selected.HH)
names(plot.data) <- c("H11", "H21", "H12", "H22")
plot.data.melted <- melt(plot.data)

ggplot(plot.data.melted, aes(x=value, color=variable)) +
  #geom_histogram(aes(y=..density..), fill="white", bins = 25, alpha=0.5, position="identity") +
  geom_density(alpha=.2)


parameter.selected.LH <- parameter.selected[which(data_simulation[,1] < low_high_boundary),]
data_simulation.LH <- data_simulation[which(data_simulation[,1] < low_high_boundary),]

parameter.selected.LH <- parameter.selected.LH[which(data_simulation.LH[,2] > low_high_boundary),]
data_simulation.LH <- data_simulation.LH[which(data_simulation.LH[,2] > low_high_boundary),]

models.LH <- rownames(parameter.selected.LH)

plot.data <- as.data.frame(parameter.selected.LH)
names(plot.data) <- c("LH11", "LH21", "LH12", "LH22")
plot.data.melted <- melt(plot.data)

ggplot(plot.data.melted, aes(x=value, color=variable)) +
  #geom_histogram(aes(y=..density..), fill="white", bins = 25, alpha=0.5, position="identity") +
  geom_density(alpha=.2)

parameter.selected.HL <- parameter.selected[which(data_simulation[,1] > low_high_boundary),]
data_simulation.HL <- data_simulation[which(data_simulation[,1] > low_high_boundary),]

parameter.selected.HL <- parameter.selected.HL[which(data_simulation.HL[,2] < low_high_boundary),]
data_simulation.HL <- data_simulation.HL[which(data_simulation.HL[,2] < low_high_boundary),]

models.HL <- rownames(parameter.selected.HL)

plot.data <- as.data.frame(parameter.selected.HL)
names(plot.data) <- c("HL11", "HL21", "HL12", "HL22")
plot.data.melted <- melt(plot.data)

ggplot(plot.data.melted, aes(x=value, color=variable)) +
  #geom_histogram(aes(y=..density..), fill="white", bins = 25, alpha=0.5, position="identity") +
  geom_density(alpha=.2)


parameter.selected <- parameters[,c(5,6,7,8)]
row.names(parameter.selected) <- seq(1:1000000)

parameter.selected <- parameter.selected[models.LH,]
plot.data <- as.data.frame(parameter.selected)
names(plot.data) <- c("L11", "L21", "L12", "L22")
plot.data.melted <- melt(plot.data)

ggplot(plot.data.melted, aes(x=value, color=variable)) +
  #geom_histogram(aes(y=..density..), fill="white", bins = 25, alpha=0.5, position="identity") +
  geom_density(alpha=.2)



#######################
# Keep only the low low state models
test <- parameters[models.LL,]
test2 <- parameters[models.HH,]
test <- rbind(test,test2)
#test[,5] <-  min(parameters[,5]) # 1.2 #
#test[,15] <- 2
#test[,14] <- 2
#test[,8] <-  max(parameters[,5]) # 1.2 #

write.table(test, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/inputs/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_0_parameters_input.txt", row.names = F , col.names = F, sep = "\t")

ic.file <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_RK_g2_ABact_g2_IC.txt"
initial.conds <- read.table(ic.file, header = F)
initial.conds <- initial.conds[models.LL,]


write.table(initial.conds, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/inputs/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_IC_input.txt", row.names = F , col.names = F, sep = "\t")

test <- parameters[models.LL,]
ic.file <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_RK_g2_ABact_g2_IC.txt"
initial.conds <- read.table(ic.file, header = F)
initial.conds <- initial.conds[models.LL,]


#######################
# Analyze processed data


output_file_base <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_1_output.txt"
data_simulation_new <- read.table(output_file_base, header = F)
#data_simulation <- data_simulation[1:50000,]
data_simulation_new <- log2(data_simulation_new)
data_simulation_new <- sweep(data_simulation_new,2,col.means,FUN = "-")
data_simulation_new <- sweep(data_simulation_new,2,col.sds,FUN = "/")


bin_count=40

plot_data <- data.frame(x=data_simulation_new[,1],y=data_simulation_new[,2])
colnames(plot_data) <- c("x", "y")
h1 <- hist(plot_data$x, breaks=bin_count, plot=F)
h2 <- hist(plot_data$y, breaks=bin_count, plot=F)
top <- max(h1$counts, h2$counts)
k <- kde2d(plot_data$x, plot_data$y, n=bin_count)

# margins
oldpar <- par()
par(mar=c(3,3,1,1))
layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
image(k, col=plot_color) #plot the image
par(mar=c(0,2,1,0))
barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
par(mar=c(2,0,0.5,1))
barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)





low_high_boundary = 0.250

parameter.file <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_1_parameters.txt"
parameters <- read.table(parameter.file, header = F)

parameters.selected <- parameters
#parameter.selected <- parameters[,c(13,14,15,16)]
row.names(data_simulation_new) <- seq(1:8513900)
test <- floor((as.numeric(rownames(data_simulation_new))-1)/100)+1
#  unlist(lapply(seq(1:85139), function(x)  rep(x,100)))


data_simulation.LL <- data_simulation_new[which(data_simulation_new[,1] < low_high_boundary),]
data_simulation.LL <- data_simulation.LL[which(data_simulation.LL[,2] < low_high_boundary),]
models.LL <- unique(floor((as.numeric(rownames(data_simulation.LL))-1)/100)+1)



data_simulation.HH <- data_simulation_new[which(data_simulation_new[,1] > low_high_boundary),]
data_simulation.HH <- data_simulation.HH[which(data_simulation.HH[,2] > low_high_boundary),]
models.HH <- unique(floor((as.numeric(rownames(data_simulation.HH))-1)/100)+1)


data_simulation.LH <- data_simulation_new[which(data_simulation_new[,1] < low_high_boundary),]
data_simulation.LH <- data_simulation.LH[which(data_simulation.LH[,2] > low_high_boundary),]
models.LH <- unique(floor((as.numeric(rownames(data_simulation.LH))-1)/100)+1)


data_simulation.HL <- data_simulation_new[which(data_simulation_new[,1] > low_high_boundary),]
data_simulation.HL <- data_simulation.HL[which(data_simulation.HL[,2] < low_high_boundary),]
models.HL <- unique(floor((as.numeric(rownames(data_simulation.LH))-1)/100)+1)

models.quad <- intersect(models.LL, models.LH)
models.quad <- intersect(models.HL, models.quad)
models.quad <- intersect(models.HH, models.quad)

parameter.selected <- parameters[models.quad,]

write.table(parameter.selected, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/inputs/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_0_parameters_input.txt", row.names = F , col.names = F, sep = "\t")

ic.file <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_RK_g2_ABact_g2_IC.txt"
initial.conds <- read.table(ic.file, header = F)
initial.conds <- initial.conds[models.quad,]


write.table(initial.conds, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/inputs/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_IC_input.txt", row.names = F , col.names = F, sep = "\t")


models <- seq(1:8513900)
selected.models <- rep(1,100*length(models.quad))
for(i in 1:length(models.quad))
{
  for(j in 1:100)
  selected.models[(i-1)*100+j] = ((models.quad[i]-1)*100 + j)
}
#test <- sapply(selected.models,function(x)(x=2))
test <- data_simulation_new[selected.models,]

bin_count = 40
plot_data <- data.frame(x=test[,1],y=test[,2])
colnames(plot_data) <- c("x", "y")
h1 <- hist(plot_data$x, breaks=bin_count, plot=F)
h2 <- hist(plot_data$y, breaks=bin_count, plot=F)
top <- max(h1$counts, h2$counts)
k <- kde2d(plot_data$x, plot_data$y, n=bin_count)

# margins
oldpar <- par()
par(mar=c(3,3,1,1))
layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
image(k, col=plot_color) #plot the image
par(mar=c(0,2,1,0))
barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
par(mar=c(2,0,0.5,1))
barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)


############################################################################################################################################

#############################################################################################################################################
#######################
# Analyze LLHH data data


output_file_base <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_1_output.txt"
data_simulation_new <- read.table(output_file_base, header = F)
#data_simulation <- data_simulation[1:50000,]
data_simulation_new <- log2(data_simulation_new)
data_simulation_new <- sweep(data_simulation_new,2,col.means,FUN = "-")
data_simulation_new <- sweep(data_simulation_new,2,col.sds,FUN = "/")



low_high_boundary = 0.250

parameter.file <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_1_parameters.txt"

parameters <- read.table(parameter.file, header = F)

ic.file <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_1_IC.txt"
ic.LLHH <-  read.table(ic.file, header = F)

parameters.selected <- parameters
#parameter.selected <- parameters[,c(13,14,15,16)]
row.names(data_simulation_new) <- seq(1:8513900)
#test <- floor((as.numeric(rownames(data_simulation_new))-1)/100)+1
#  unlist(lapply(seq(1:85139), function(x)  rep(x,100)))


data_simulation.LL <- data_simulation_new[which(data_simulation_new[,1] < low_high_boundary),]
data_simulation.LL <- data_simulation.LL[which(data_simulation.LL[,2] < low_high_boundary),]
models.LL <- unique(floor((as.numeric(rownames(data_simulation.LL))-1)/100)+1)



data_simulation.HH <- data_simulation_new[which(data_simulation_new[,1] > low_high_boundary),]
data_simulation.HH <- data_simulation.HH[which(data_simulation.HH[,2] > low_high_boundary),]
models.HH <- unique(floor((as.numeric(rownames(data_simulation.HH))-1)/100)+1)


data_simulation.LH <- data_simulation_new[which(data_simulation_new[,1] < low_high_boundary),]
data_simulation.LH <- data_simulation.LH[which(data_simulation.LH[,2] > low_high_boundary),]
models.LH <- unique(floor((as.numeric(rownames(data_simulation.LH))-1)/100)+1)


data_simulation.HL <- data_simulation_new[which(data_simulation_new[,1] > low_high_boundary),]
data_simulation.HL <- data_simulation.HL[which(data_simulation.HL[,2] < low_high_boundary),]
models.HL <- unique(floor((as.numeric(rownames(data_simulation.LH))-1)/100)+1)

#######################
# Keep quad models

models.quad <- intersect(models.LL, models.LH)
models.quad <- intersect(models.HL, models.quad)
models.quad <- intersect(models.HH, models.quad)

parameter.selected <- parameters[models.quad,]

write.table(parameter.selected, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/quad_parameters.txt", row.names = F , col.names = F, sep = "\t")


numberIC <- 100
quad.models.all.ic <- rep(models.quad -1,each = numberIC)
for(i in 1:length(quad.models.all.ic)){
  quad.models.all.ic[i] <- 100*quad.models.all.ic[i] + 1 + (i-1) %% numberIC
}

data_simulation.quad <- data_simulation_new[quad.models.all.ic,]

write.table(data_simulation.quad, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/quad_output_logNormal.txt", row.names = F , col.names = F, sep = "\t")

data_simulation.quad.LL <- data_simulation.quad[which(data_simulation.quad[,1] < low_high_boundary),]
data_simulation.quad.LL <- data_simulation.quad.LL[which(data_simulation.quad.LL[,2] < low_high_boundary),]
data_simulation.quad.LL.ic <- rownames(data_simulation.quad.LL)

data_simulation.quad.LL.ic <- ic.LLHH[data_simulation.quad.LL.ic,]
write.table(data_simulation.quad.LL.ic, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/quad_LL_IC.txt", row.names = F , col.names = F, sep = "\t")


models.multiple.quad.LL.ic <- (floor((as.numeric(rownames(data_simulation.quad.LL.ic))-1)/100)+1)
parameter.multiple.quad.LL.ic <- parameters[models.multiple.quad.LL.ic,]

write.table(parameter.multiple.quad.LL.ic, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/quad_parameters_LL_IC_multiple.txt", row.names = F , col.names = F, sep = "\t")

write.table(data_simulation.quad.LL.ic, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/inputs/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_0_IC_input.txt", row.names = F , col.names = F, sep = "\t")

tmp.parameters <- parameter.multiple.quad.LL.ic
#reduction_fator <- 0.97
#tmp.parameters[,5] <- (reduction_fator/1.5)*tmp.parameters[,5]
#tmp.parameters[,8] <- (reduction_fator)*tmp.parameters[,8]


#tmp.parameters[,7] <- 5*tmp.parameters[,7]
tmp.parameters[,15] <- max(1,0.05*tmp.parameters[,15])
#tmp.parameters[,13] <- max(1,(reduction_fator/1.5)*tmp.parameters[,13])
#tmp.parameters[,16] <- max(1,(reduction_fator)*tmp.parameters[,16])


write.table(tmp.parameters, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/inputs/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_0_parameters_input.txt", row.names = F , col.names = F, sep = "\t")

setwd("/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic")

output_file <- sRACIPE_simulate_GRN(topology_file = topology_file, MAX_NOISE = 0, NOISE_LEVELS = 1, NUM_MODELS = 10486, INITIAL_CONDITIONS = 1, ANNEAL=F, NOISE_SCALING_FACTOR=0.6, SIM_TIME = 50, STEP_SIZE = 0.01,PARAMETERS_FILE = T,  READ_IC = T)

output_file_base <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_EM_g2_ABact_g2_Annealing_1_pmFile_1_output.txt"

data_simulation_quad_LL <- read.table(output_file_base, header = F)

data_simulation_test <- data_simulation_quad_LL[,61:62]
data_simulation_test <- log2(data_simulation_test)
data_simulation_test <- sweep(data_simulation_test,2,col.means,FUN = "-")
data_simulation_test <- sweep(data_simulation_test,2,col.sds,FUN = "/")

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


########
# Keep only the models without repetition
#rownames(tmp.parameters) <- seq(1:dim(tmp.parameters)[1])
#test <- rownames(unique(tmp.parameters[,1:2]))
#tmp.parameters <- tmp.parameters[test,]


output_file <- sRACIPE_simulate_GRN(topology_file = topology_file, MAX_NOISE = 50, NOISE_LEVELS = 50, NUM_MODELS = 350, INITIAL_CONDITIONS = 100, ANNEAL=T, NOISE_SCALING_FACTOR=0.8, SIM_TIME = 500, STEP_SIZE = 0.01,PARAMETERS_FILE = T,  READ_IC = F)

output_file <- sRACIPE_simulate_GRN(topology_file = topology_file, MAX_NOISE = 50, NOISE_LEVELS = 50, NUM_MODELS = 350, INITIAL_CONDITIONS = 100, ANNEAL=T, NOISE_SCALING_FACTOR=0.8, SIM_TIME = 500, STEP_SIZE = 0.01,PARAMETERS_FILE = T,  READ_IC = F)


######################################################################################
))
stop()


output_file_base <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_EM_g2_ABact_g2_Annealing_1_pmFile_1_output.txt"

data_simulation_quad_LL <- read.table(output_file_base, header = F)
data_simulation_test <- data_simulation_quad_LL[,99:100]
data_simulation_test <- log2(data_simulation_test)
data_simulation_test <- sweep(data_simulation_test,2,col.means,FUN = "-")
data_simulation_test <- sweep(data_simulation_test,2,col.sds,FUN = "/")

test <- data_simulation_test


bin_count = 40
plot_data <- data.frame(x=test[,1],y=test[,2])
colnames(plot_data) <- c("x", "y")

ggplot(plot_data, aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #geom_contour((aes(x=x, y=y,z = ..density..))) +
  scale_fill_distiller(palette= "Spectral", direction = -1)  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlim(-2.5,2) +
  ylim(-2.5,2) +
  labs(x = 'Normalized A expression('~log[2]~')', y = 'Normalized B expression('~log[2]~')') +
  theme(
    legend.position='none',
    text = element_text(size=20)
  )




h1 <- hist(plot_data$x, breaks=bin_count, plot=F)
h2 <- hist(plot_data$y, breaks=bin_count, plot=F)
top <- max(h1$counts, h2$counts)
k <- kde2d(plot_data$x, plot_data$y, n=bin_count)
#image(k, col=plot_color, xlim = c(-2.5,2), ylim = c(-2.5,2)) #plot the image

# margins
oldpar <- par()
par(mar=c(3,3,1,1))
layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
image(k, col=plot_color) #plot the image
par(mar=c(0,2,1,0))
barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
par(mar=c(2,0,0.5,1))
barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)



############################################################################################################################################

############################################################################################################################################
# Keep only the LL models
data_simulation_LL <- data_simulation_new[1:6044700,]

data_simulation.LL <- data_simulation_LL[which(data_simulation_LL[,1] < low_high_boundary),]
data_simulation.LL <- data_simulation.LL[which(data_simulation.LL[,2] < low_high_boundary),]
models.LL <- unique(floor((as.numeric(rownames(data_simulation.LL))-1)/100)+1)



data_simulation.HH <- data_simulation_LL[which(data_simulation_LL[,1] > low_high_boundary),]
data_simulation.HH <- data_simulation.HH[which(data_simulation.HH[,2] > low_high_boundary),]
models.HH <- unique(floor((as.numeric(rownames(data_simulation.HH))-1)/100)+1)


data_simulation.LH <- data_simulation_LL[which(data_simulation_LL[,1] < low_high_boundary),]
data_simulation.LH <- data_simulation.LH[which(data_simulation.LH[,2] > low_high_boundary),]
models.LH <- unique(floor((as.numeric(rownames(data_simulation.LH))-1)/100)+1)


data_simulation.HL <- data_simulation_LL[which(data_simulation_LL[,1] > low_high_boundary),]
data_simulation.HL <- data_simulation.HL[which(data_simulation.HL[,2] < low_high_boundary),]
models.HL <- unique(floor((as.numeric(rownames(data_simulation.LH))-1)/100)+1)

models.quad <- intersect(models.LL, models.LH)
models.quad <- intersect(models.HL, models.quad)
models.quad <- intersect(models.HH, models.quad)

selected.models <- rep(1,100*length(models.quad))
for(i in 1:length(models.quad))
{
  for(j in 1:100)
    selected.models[(i-1)*100+j] = ((models.quad[i]-1)*100 + j)
}

#test <- sapply(selected.models,function(x)(x=2))

#test <- sapply(selected.models,function(x)(x=2))
test <- data_simulation_new[selected.models,]

bin_count = 40
plot_data <- data.frame(x=test[,1],y=test[,2])
colnames(plot_data) <- c("x", "y")
h1 <- hist(plot_data$x, breaks=bin_count, plot=F)
h2 <- hist(plot_data$y, breaks=bin_count, plot=F)
top <- max(h1$counts, h2$counts)
k <- kde2d(plot_data$x, plot_data$y, n=bin_count)

# margins
oldpar <- par()
par(mar=c(3,3,1,1))
layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
image(k, col=plot_color) #plot the image
par(mar=c(0,2,1,0))
barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
par(mar=c(2,0,0.5,1))
barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)

setwd("/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic")

output_file <- sRACIPE_simulate_GRN(topology_file = topology_file, PARAMETERS_FILE = 1, MAX_NOISE = 30, NOISE_LEVELS = 50, NUM_MODELS = 100, INITIAL_CONDITIONS = 10, ANNEAL=T, NOISE_SCALING_FACTOR=0.9, SIM_TIME = 100, STEP_SIZE = 0.05)

output_file_base <- "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/results/sRACIPE_EM_g2_ABact_g2_Annealing_1_pmFile_1_output.txt"

data_simulation_LL_anneal <- read.table(output_file_base, header = F)
data_simulation_test <- data_simulation_LL_anneal[,99:100]
data_simulation_test <- log2(data_simulation_test)
data_simulation_test <- sweep(data_simulation_test,2,col.means,FUN = "-")
data_simulation_test <- sweep(data_simulation_test,2,col.sds,FUN = "/")


bin_count=40

plot_data <- data.frame(x=data_simulation_test[,1],y=data_simulation_test[,2])
colnames(plot_data) <- c("x", "y")
h1 <- hist(plot_data$x, breaks=bin_count, plot=F)
h2 <- hist(plot_data$y, breaks=bin_count, plot=F)
top <- max(h1$counts, h2$counts)
k <- kde2d(plot_data$x, plot_data$y, n=bin_count)

# margins
oldpar <- par()
par(mar=c(3,3,1,1))
layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
image(k, col=plot_color) #plot the image
par(mar=c(0,2,1,0))
barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
par(mar=c(2,0,0.5,1))
barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)

########################################

test <- floor((as.numeric(rownames(data_simulation_new))-1)/100)+1

plot.data <- as.data.frame(parameter.selected[,13:16])
plot.data.melted <- melt(plot.data)
ggplot(plot.data.melted, aes(x=value, color=variable)) +
  geom_density(alpha=.2)


parameter.selected.HH <- parameter.selected[which(data_simulation[,1] > low_high_boundary),]
data_simulation.HH <- data_simulation[which(data_simulation[,1] > low_high_boundary),]

parameter.selected.HH <- parameter.selected.HH[which(data_simulation.HH[,2] > low_high_boundary),]
data_simulation.HH <- data_simulation.HH[which(data_simulation.HH[,2] > low_high_boundary),]

models.HH <- rownames(parameter.selected.HH)

plot.data <- as.data.frame(parameter.selected.HH)
names(plot.data) <- c("H11", "H21", "H12", "H22")
plot.data.melted <- melt(plot.data)

ggplot(plot.data.melted, aes(x=value, color=variable)) +
  #geom_histogram(aes(y=..density..), fill="white", bins = 25, alpha=0.5, position="identity") +
  geom_density(alpha=.2)


parameter.selected.LH <- parameter.selected[which(data_simulation[,1] < low_high_boundary),]
data_simulation.LH <- data_simulation[which(data_simulation[,1] < low_high_boundary),]

parameter.selected.LH <- parameter.selected.LH[which(data_simulation.LH[,2] > low_high_boundary),]
data_simulation.LH <- data_simulation.LH[which(data_simulation.LH[,2] > low_high_boundary),]

models.LH <- rownames(parameter.selected.LH)

plot.data <- as.data.frame(parameter.selected.LH)
names(plot.data) <- c("LH11", "LH21", "LH12", "LH22")
plot.data.melted <- melt(plot.data)

ggplot(plot.data.melted, aes(x=value, color=variable)) +
  #geom_histogram(aes(y=..density..), fill="white", bins = 25, alpha=0.5, position="identity") +
  geom_density(alpha=.2)

parameter.selected.HL <- parameter.selected[which(data_simulation[,1] > low_high_boundary),]
data_simulation.HL <- data_simulation[which(data_simulation[,1] > low_high_boundary),]

parameter.selected.HL <- parameter.selected.HL[which(data_simulation.HL[,2] < low_high_boundary),]
data_simulation.HL <- data_simulation.HL[which(data_simulation.HL[,2] < low_high_boundary),]

models.HL <- rownames(parameter.selected.HL)

plot.data <- as.data.frame(parameter.selected.HL)
names(plot.data) <- c("HL11", "HL21", "HL12", "HL22")
plot.data.melted <- melt(plot.data)

ggplot(plot.data.melted, aes(x=value, color=variable)) +
  #geom_histogram(aes(y=..density..), fill="white", bins = 25, alpha=0.5, position="identity") +
  geom_density(alpha=.2)


parameter.selected <- parameters[,c(5,6,7,8)]
row.names(parameter.selected) <- seq(1:1000000)

parameter.selected <- parameter.selected[models.LH,]
plot.data <- as.data.frame(parameter.selected)
names(plot.data) <- c("L11", "L21", "L12", "L22")
plot.data.melted <- melt(plot.data)

ggplot(plot.data.melted, aes(x=value, color=variable)) +
  #geom_histogram(aes(y=..density..), fill="white", bins = 25, alpha=0.5, position="identity") +
  geom_density(alpha=.2)


write.table(parameter.selected, file = "/Users/koharv/Documents/Work/sRACIPE_project/eLife_synthetic/inputs/sRACIPE_EM_g2_ABact_g2_Annealing_0_pmFile_0_parameters_input.txt", row.names = F , col.names = F, sep = "\t")




##########################################################################################################################################

,xlim=c(-2.5, 1.50),ylim=c(-2.5, 1.50)

##########################################################################################################################################
ggplot(df, aes(x=weight)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")

ggplot(plot.data.melted, aes(x=plot.data))
       + geom_density()

hist(parameter.selected.LL[,4], breaks = 40)




#rm(list = ls())

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


working_directory <- getwd()
config_file = "/Users/koharv/Documents/Work/sRACIPE_dev/inputs/sRACIPE.cfg"
configuration <- sRACIPE_load_configuration(config_file = config_file)

working_directory <- getwd()
output_file_base <- "/Users/koharv/Documents/Work/sRACIPE_dev/results/sRACIPE_RK_g2_ABact_g2_output.txt"
name_genes <- read.table(paste(working_directory,"/results/gene_interaction_topology_",topology$filename,".txt",sep=""), header = T, stringsAsFactors = F)
data_simulation <- read.table(output_file_base, header = F)
#data_simulation <- data_simulation[1:50000,]
data_simulation <- log2(data_simulation)
name_genes <- t(as.matrix(name_genes))
col.means <- colMeans(data_simulation)
col.sds <- apply(data_simulation, 2, sd)
data_simulation <- sweep(data_simulation,2,col.means,FUN = "-")
data_simulation <- sweep(data_simulation,2,col.sds,FUN = "/")

#parameter.selected <- parameters[,c(13,15)]
#parameter.selected <- rep(parameter.selected, each = 100)
#data.selected <- data.simulation[which(parameter.selected > 90),]








topology_file <- "/Users/koharv/Documents/Work/sRACIPE_dev/inputs/g2_ABactLL.tpo"
topology <- sRACIPE_load_topology(topology_file = topology_file)
output_file <- sRACIPEv03::sRACIPE_RK_deterministic_thMod(topology_file = topology_file, NUM_MODELS = 2000, FCH_MAX = 1.1, FCH_MIN = 1)
data_simulation1 <- read.table(output_file, header = F)
data_simulation1 <- log2(data_simulation1)
#col.means <- colMeans(data_simulation)
#col.sds <- apply(data_simulation, 2, sd)

data_simulation1 <- sweep(data_simulation1,2,col.means,FUN = "-")
data_simulation1 <- sweep(data_simulation1,2,col.sds,FUN = "/")





library(reshape2)
library(ggplot2)
library(Rtsne)
library(gplots)
library(MASS)
library(RColorBrewer)
#library(htmlwidgets)
#library(d3heatmap)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot_color <- rf(32)


output_file <- "/Users/koharv/Documents/Work/sRACIPE_dev/results/sRACIPE_RK_test_g2_output.txt"
topology_file <- "/Users/koharv/Documents/Work/sRACIPE_dev/inputs/test.tpo"
topology <- sRACIPE_load_topology(topology_file = topology_file)

initial.conditions <- 100
num.models <- 10000
data.simulation <- read.table(output_file, header = F)
data.simulation <- load_data(data.simulation, topology)


bin.count <- 40
plot.filename <- "plots"


#pca1 = prcomp(data.simulation, scale. = FALSE)
#pca.data <- data.frame(x=pca1$x[,1],y=pca1$x[,2])
pca.data <- data.frame(x=data.simulation[,1],y=data.simulation[,2])


pdf(paste("results/",plot.filename,"_deterministic.pdf",sep = ""))

density_plot(pca.data,bin.count,plot_color)
dev.off()

parameter.file <- "/Users/koharv/Documents/Work/sRACIPE_dev/results/sRACIPE_RK_test_g2_parameters.txt"
working_directory <- getwd()
parameters <- read.table(parameter.file, header = F)

#parameter.selected <- parameters[,c(13,15)]
#parameter.selected <- rep(parameter.selected, each = 100)
#data.selected <- data.simulation[which(parameter.selected > 90),]


parameter.selected <- parameters[,c(13,14,15,16)]
row.names(parameter.selected) <- seq(1:num.models)
parameter.selected <- parameter.selected[which(parameter.selected[,1] < 20),]
parameter.selected <- parameter.selected[which(parameter.selected[,2] < 20),]
parameter.selected <- parameter.selected[which(parameter.selected[,3] < 20),]
parameter.selected <- parameter.selected[which(parameter.selected[,4] < 20),]
parameter.selected <- as.integer(row.names(parameter.selected))
tmp <- rep(seq(-(initial.conditions-1),0),length(parameter.selected))
parameter.selected <- rep(parameter.selected, each = initial.conditions)
parameter.selected <- 100*parameter.selected + tmp

data.selected <- data.simulation[parameter.selected,]

pdf(paste("results/",plot.filename,"_lowA_lowB_deterministic.pdf",sep = ""))
pca.data <- data.frame(x=data.selected[,1],y=data.selected[,2])
density_plot(pca.data,bin.count,plot_color)
dev.off()

data.ic.col <- matrix(data.simulation,nrow = num.models, byrow = F)

data_simulation <- scale(data_simulation)
name_models <- seq(1:nrow(data_simulation))
row.names(data_simulation) <- name_models
colnames(data_simulation) <- name_genes

data.test <- unique(floor((which((data.simulation[,1] > .25) & (data.simulation[,2] > .25))-1)/100))

data.test2 <- unique(floor((which((data.simulation[,1] < .25) & (data.simulation[,2] < .25))-1)/100))

parameter.selected <- parameters[intersect(data.test, data.test2),]

parameter.selected1 <- parameter.selected[which((data.simulation[,1] < 0.2) && (data.simulation[,2] < 0.2))]
row.names(parameter.selected) <- seq(1:num.models)
parameter.selected <- parameter.selected[which(parameter.selected[,1] < 20),]
parameter.selected <- parameter.selected[which(parameter.selected[,2] < 20),]
parameter.selected <- parameter.selected[which(parameter.selected[,3] < 20),]
parameter.selected <- parameter.selected[which(parameter.selected[,4] < 20),]
parameter.selected <- as.integer(row.names(parameter.selected))
tmp <- rep(seq(-(initial.conditions-1),0),length(parameter.selected))
parameter.selected <- rep(parameter.selected, each = initial.conditions)
parameter.selected <- 100*parameter.selected + tmp

data.selected <- data.simulation[parameter.selected,]

dataA <- matrix(data.simulation[,1],nrow = num.models, ncol = initial.conditions, byrow = T)
test1 <- 10*dataA
test1 <- round(test1)
num.state.A <- (apply(test1,1,unique))
test2 <- (lapply(num.state.A,length))
test3 <-

  dim(num.state.A)
num.state.A[num.state.A<0.01] <- 0
num.state.A <- t(apply(num.state.A,1,unique))





topology_file <- "/Users/koharv/Documents/Work/sRACIPE_dev/inputs/g2_ABactLL.tpo"
topology <- sRACIPE_load_topology(topology_file = topology_file)
output_file <- sRACIPEv03::sRACIPE_RK_deterministic_thMod(topology_file = topology_file, NUM_MODELS = 2000, FCH_MAX = 1.1, FCH_MIN = 1)
data_simulation1 <- read.table(output_file, header = F)
data_simulation1 <- log2(data_simulation1)
#col.means <- colMeans(data_simulation)
#col.sds <- apply(data_simulation, 2, sd)

data_simulation1 <- sweep(data_simulation1,2,col.means,FUN = "-")
data_simulation1 <- sweep(data_simulation1,2,col.sds,FUN = "/")

bin_count=40


plot_data <- data.frame(x=data_simulation1[,1],y=data_simulation1[,2])
colnames(plot_data) <- c("x", "y")
h1 <- hist(plot_data$x, breaks=bin_count, plot=F)
h2 <- hist(plot_data$y, breaks=bin_count, plot=F)
top <- max(h1$counts, h2$counts)
k <- kde2d(plot_data$x, plot_data$y, n=bin_count)

# margins
oldpar <- par()
par(mar=c(3,3,1,1))
layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
image(k, col=plot_color) #plot the image
par(mar=c(0,2,1,0))
barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
par(mar=c(2,0,0.5,1))
barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)


