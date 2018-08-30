library(GEOquery)
library(Biobase)
library(scran)
library(SingleCellExperiment)
library(scater)
#library(biomart)
#gse <- getGEO("GSE110357", GSEMatrix = TRUE)
#show(gse)
#test <- pData(gse)
#dim(pData(gse[[1]]))
#head(Meta(gse))
#names(GSMList(gse))
#show(pData(phenoData(gse[[1]]))[1:5,c(1,6,8)])

#show(pData(phenoData(gse[[1]]))[1:5,c(1,6,8)])

#test <- exprs(gse[[1]])

#gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
#head(gsmplatforms)

#eset <- pData(gse[[1]])

#gds <- getGEO("GSE110357")
#eset <- GDS2eSet(gse)


data_nature <- read.csv("data/GSE110357_htseq_counts_all_v1.csv",header = T)
rownames(data_nature) <- data_nature[,1]
data_nature <- data_nature[,2:384]
dim(data_nature)
test <- rowSums(data_nature)
data_nature <- data_nature[which(test>20),]

sce <- SingleCellExperiment(list(counts=as.matrix((data_nature))))
sce <- computeSumFactors(sce)
#summary(sizeFactors(sce))
sce <- normalize(sce)
umi <-
  plotPCA(
    sce )

#lc_sce <-logcounts(sce)

pca_data <- data.frame(x=umi$data$X,y=umi$data$Y)
ggplot(pca_data, aes(x=x, y=y) ) +
  #geom_point(aes(color=z))
  geom_point( alpha=0.5)

gene_names <- c("Zeb1", "Cdh2", "Cdh11", "Krt5", "Col3a1", "Vim", "Cdh1", "Epcam", "Snai1" ,"Smad2", "Esrp1", "Twist1", "Trp63", "Klf5", "Cebpa" ,"Grhl2")

gene <- list("ENSMUSG00000024238","ENSMUSG00000024304","ENSMUSG00000031673","ENSMUSG00000061527","ENSMUSG00000026043","ENSMUSG00000026728","ENSMUSG00000000303","ENSMUSG00000045394","ENSMUSG00000042821","ENSMUSG00000024563","ENSMUSG00000040728","ENSMUSG00000035799","ENSMUSG00000022510","ENSMUSG00000005148","ENSMUSG00000034957","ENSMUSG00000022286")
getData <- function(data_nature,gene){return(data_nature[grep(gene, rownames(data_nature)),])}
geneData <- lapply(gene,function(X) getData(data_nature,X))
geneData <- matrix(unlist(geneData), ncol = 383, byrow = T)
geneData <- geneData + 1
geneData <- log2(geneData)
rownames(geneData) <- gene_names
heatmap((geneData), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), cexRow = 2)
#tmp <- scale(geneData)
#heatmap((tmp), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), cexRow = 2)
TFstoKeep <- c ("Krt5", "Klf5", "Cdh1","Esrp1", "Epcam", "Trp63", "Cebpa" ,"Grhl2", "Twist1","Snai1" , "Zeb1", "Cdh2","Smad2",   "Cdh11", "Vim",  "Col3a1")

pca_expData = prcomp(t(geneData),center = T, scale. =T)
pca_data <- data.frame(x=pca_expData$x[,1],y=pca_expData$x[,2])

summary(pca_expData)

ggplot(pca_data, aes(x=x, y=y) ) +
  #geom_point(aes(color=z))
  geom_point( alpha=0.5)



p = list()
for(i in 1:4){#length(TFstoKeep)) {
  title = TFstoKeep[i]
  pca_data$z = geneData[title,]
  p[[i]] = plotScatter(pca_data, title)
}
grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]],  ncol=4)

do.call(grid.arrange,p)

#lay <- c(1,2,3,4)



library(gridExtra)
plotScatter <- function(pca_data, title){
  (ggplot(pca_data, aes(x=x, y=y), aspect ) +
     #geom_point(aes(color=z))
     geom_point(aes(color=z), alpha=0.5) +
     scale_color_gradient(low = "darkblue", high = "orange", name = title) +
     xlim(-8,5) +
      ylim(-5,5) +
     #labs(x = 'PC1', y = 'PC2') +
     #guides(fill=guide_legend(title=title)) +
     theme(
       legend.position='none',
       # legend.title = element_text(title), #element_text(colour="blue", size=10, face="bold"),
       text = element_text(size=15),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       axis.title.x=element_blank(),
       #axis.text.x=element_blank(),
       axis.title.y=element_blank(),
       #axis.text.y=element_blank(),
       #axis.ticks.x=element_blank(),
       axis.line = element_line(colour = "black"),
       aspect.ratio = 1
     ))
}

ggplot(pca_data, aes(x=x, y=y) ) +
    #geom_point(aes(color=z))
    geom_point(aes(color=z), alpha=0.5) +
    scale_color_gradient(low = "darkblue", high = "orange", name = title) +
    #xlim(-6,4.5) +
    # ylim(-2.5,2) +
    #labs(x = 'PC1', y = 'PC2') +
    #guides(fill=guide_legend(title=title)) +
    theme(
      legend.position='top',
      # legend.title = element_text(title), #element_text(colour="blue", size=10, face="bold"),
      text = element_text(size=15),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      axis.title.y=element_blank(),
      #axis.text.y=element_blank(),
      #axis.ticks.x=element_blank(),
      axis.line = element_line(colour = "black")
    )


##############################################

pca_data <- scale((tmp)) %*% pca_all_exp$rotation
pca_data <- data.frame(x=pca_all_exp[,1],y=pca_all_exp[,2])



data_plot <- log2(data_nature+1)


geneNames <- rownames(test)
library("AnnotationDbi")
library("org.Mm.eg.db")
res<-results(dds,alpha=.05, contrast=c("Type", "Disease", "Control"))
res$symbol <- mapIds(org.Hs.eg.db,keys=row.names(res),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
#gene_names <- c("Zeb1" "Cdh2" "Cdh11" "Krt5" "Col3a1" "Vim" "Cdh1" "Epcam" "Snai1" "Smad2" "Esrp1" "Twist1" "Trp63" "Klf5" "Cebpa" "Grhl2")

test2 <- test


plotPCA(
  sce,ntop=500)

tmp <- scale(geneData)
pca1 = prcomp(t(tmp), scale. =F)

pca_data <- data.frame(x=pca1$x[,1],y=pca1$x[,2])

rownames(geneData) <- name_genes.modified
rownames(tmp) <- name_genes.modified

heatmap((tmp), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2),cexRow = 2)

Epcam=tmp["Epcam",]
#names(z) <- tmpGene
ggplot(pca_data, aes(x=x, y=y) ) +
  #geom_point(aes(color=z))
 geom_point(aes(color=pca_data$z), alpha = 0.75, size=4) +
  scale_color_gradient(low = "darkblue", high = "orange",breaks = c(-1, 0, 1), labels = c(-1, 0, 1)) +
 # scale_color_gradient(low = "darkblue", high = "orange",  limits=c(-1.,1.)) +
  #scale_fill_distiller(palette= "Spectral", direction = -1) +

  #stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #geom_contour((aes(x=x, y=y,z = ..density..))) +
  #scale_fill_distiller(palette= "Spectral", direction = -1)  +
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  #xlim(-10,5) +
  #ylim(-3,3) +
  labs(x = 'PC1', y = 'PC2') +
  theme(
    legend.position='top',
    #legend.title = element_blank(), #element_text(colour="blue", size=10, face="bold"),
    text = element_text(size=45),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


Klf=tmp["Klf",]
#names(z) <- tmpGene
ggplot(pca_data, aes(x=x, y=y) ) +
  #geom_point(aes(color=z))
  geom_point(aes(color=Klf), size=6) +
  scale_color_gradient(low = "darkblue", high = "orange",breaks = c(-1, 0, 1), labels = c(-1, 0, 1)) +
  # scale_color_gradient(low = "darkblue", high = "orange",  limits=c(-1.,1.)) +
  #scale_fill_distiller(palette= "Spectral", direction = -1) +

  #stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #geom_contour((aes(x=x, y=y,z = ..density..))) +
  #scale_fill_distiller(palette= "Spectral", direction = -1)  +
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  #xlim(-10,5) +
  #ylim(-3,3) +
  labs(x = 'PC1', y = 'PC2') +
  theme(
    legend.position='top',
    #legend.title = element_blank(), #element_text(colour="blue", size=10, face="bold"),
    text = element_text(size=45),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


pca_all_exp <- prcomp(t(test), scale. =T)
pca_data <- data.frame(x=pca_all_exp$x[,1],y=pca_all_exp$x[,2])

