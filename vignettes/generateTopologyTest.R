

# Generate three noede topologies by inorporating inhibitory links
# in the original 13 topologies from Uri Elon's book.
#counter <- 13
GenerateTopologies <- function(){
  l1.tpos <- seq(1,13)
  for(l1.tpo.counter in 1:length(l1.tpos))
  {
    #l1.tpo.counter <- 1 # test
    tpo.filename <- paste("g3_",l1.tpos[l1.tpo.counter],"_1.txt", sep = "")
    topology.file <- file.path(getwd(),"inputs",tpo.filename)

    f.tpo <- read.table(topology.file, header = TRUE,stringsAsFactors = FALSE)
    links <- dim(f.tpo)[1]
    links <- 5
    inhib.file <-  2^links
    new.type <- matrix(0, nrow = inhib.file, ncol = links, byrow = F )
    for(link.counter in 1: links)
    {
      new.type[,link.counter] <- rep(1:2, each = 2^(link.counter-1), len = inhib.file)
    }

    colnames(f.tpo) <- c("Source","Target","Type")
    for(inhib.file.counter in 1: inhib.file)
    {
      #f.tpo$Type <- new.type[inhib.file.counter,]
      print(new.type[inhib.file.counter,])
      tpo.filename.tmp <- paste("g3_",l1.tpos[l1.tpo.counter],"_",inhib.file.counter,".txt", sep = "")
      #write.table(f.tpo, tpo.filename.tmp, quote = FALSE, row.names = FALSE)
      #counter <- counter +1

    }

  }
}


#counter


#counter
