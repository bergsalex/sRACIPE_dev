library(GEOquery)
library(Biobase)
library(scran)
library(SingleCellExperiment)
library(scater)
gse <- getGEO("GSE110357", GSEMatrix = TRUE)
show(gse)
test <- pData(gse)
dim(pData(gse[[1]]))
head(Meta(gse))
names(GSMList(gse))
show(pData(phenoData(gse[[1]]))[1:5,c(1,6,8)])

show(pData(phenoData(gse[[1]]))[1:5,c(1,6,8)])

test <- exprs(gse[[1]])

gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
head(gsmplatforms)

eset <- pData(gse[[1]])

gds <- getGEO("GSE110357")
eset <- GDS2eSet(gse)


data_nature <- read.csv("data/GSE110357_htseq_counts_all_v1.csv",header = T)
rownames(data_nature) <- data_nature[,1]
data_nature <- data_nature[,2:384]
dim(data_nature)
test <- rowSums(data_nature)
data_nature <- data_nature[which(test>20),]

data_plot <- log2(data_nature+1)

pca1 = prcomp(t(test), scale. =F)

pca_data <- data.frame(x=pca1$x[,1],y=pca1$x[,2])


ggplot(pca_data, aes(x=x, y=y) ) +
  geom_point()
+
  #stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #geom_contour((aes(x=x, y=y,z = ..density..))) +
  scale_fill_distiller(palette= "Spectral", direction = -1)  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlim(-15,20) +
  ylim(-15,20) +
  labs(x = 'Normalized Y expression('~log[2]~')', y = 'Normalized X expression('~log[2]~')') +
  theme(
    legend.position='none',
    text = element_text(size=20)
  )





sce <- SingleCellExperiment(list(counts=as.matrix((data_nature))))
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))
sce <- normalize(sce)
test <-logcounts(sce)

umi <-
  plotPCA(
  sce, ntop = 500)



