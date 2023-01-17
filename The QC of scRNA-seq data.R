#QC and selecting cells for further analysis
library(Seurat)
args = commandArgs(trailingOnly=TRUE)

readdata <- function(dir,name){
	print(paste(dir,'/',name,sep=""))
	a1.data <- Read10X(data.dir = paste(dir,'/',name,sep=""))
	colnames(a1.data)= paste(name,colnames(a1.data),sep="")
	a1 <- CreateSeuratObject(raw.data = a1.data)
	a1@meta.data$group <- name
	return (a1)
}

calc.filter <- function(c1){
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = c1@raw.data), value = TRUE)
	percent.mito <- Matrix::colSums(c1@raw.data[mito.genes, ])/Matrix::colSums(c1@raw.data)
	c1 <- AddMetaData(object = c1, 
		metadata = percent.mito, 
		col.name = "percent.mito")
	c1 <- FilterCells(object = c1, 
		subset.names = c("nGene","percent.mito"),
		low.thresholds=c(200,-Inf),
		high.thresholds = c(4000,0.05))
	c1 <- SubsetData(object = c1,cells.use= c1@cell.names,subset.raw = T,min.genes = 10)
	return (c1)
}
calc.runpca <- function(x, jack = TRUE) {
	x <- NormalizeData(x)
	x <- FindVariableGenes(x, do.plot = FALSE)
	x <- ScaleData(x)
	x <- RunPCA(x, pc.genes = x@var.genes,
                     pcs.compute = 70, do.print = FALSE)
	return(x)
}

calc.tsne <- function(x, pcs, res=1.2) {
	x <- FindClusters(x, dims.use = 1:pcs, resolution = res,
                    print.output = FALSE, save.SNN = TRUE)
	x <- RunTSNE(x, dims.use = 1:pcs, check_duplicates = FALSE,
        		do.fast = TRUE)
	return(x)
}
calc.findpcsnum <- function(x,thres=0.001){
	pcs.num<-70
	x <- JackStraw(x, num.pc = 70, do.par = TRUE, num.cores = 5)
	pcs.factor<-levels(JackStrawPlot(x,PCs = 1:50)$data$PC.Score)
	for(i in 1:length(pcs.factor))
	{
        	print(i)
	        tt<-as.character(pcs.factor)[i]
        	y<-strsplit(tt,' ')[[1]]
	        if(as.numeric(y[2])>thres){
	                pcs.num=i
	                break
	        }
	}
	return (pcs.num)
}
calc.findpcsnum1 <- function(x,thres=0.001){
        pcs.num<-70
        x <- JackStraw(x, num.pc = 70, do.par = TRUE, num.cores = 5)
        x <- JackStrawPlot(x,PCs = 1:70)
        y <- x@dr$pca@jackstraw@overall.p.values[,2]
        print(y)
        for(i in 1:length(y))
        {
                if(y[i]>thres){
                        pcs.num=i-1
                        break
                }
        }
        return (pcs.num)
}

namelist<-read.table(args[1])
name_list<-as.vector(namelist[[1]])
indir<-"."
outdir<-"./merged"
dat<-vector("list",dim(namelist)[1])
id=1
for(name in name_list){
	print(name)
	print("reading data......")
	dat[[id]] <- readdata(indir,name)
#	dat[[id]] <- calc.filtgenes(dat[[id]])
        if(id==1)
                mdat<-dat[[id]]
        else
                mdat<-MergeSeurat(object1=mdat,object2=dat[[id]])
        id<-id+1
}

#calculate the expression percent of mitochondrial genes
mito.genes<-c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")
percent.mito <- Matrix::colSums(mdat@raw.data[mito.genes, ])/Matrix::colSums(mdat@raw.data)
mdat@meta.data$mt<-percent.mito
mdat<-FilterCells(object = mdat,
                subset.names = c("nGene","mt","nUMI"),
                low.thresholds=c(200,-Inf,1000),
                high.thresholds = c(5000,0.1,40000))
mdat<-SubsetData(mdat,cells.use=mdat@cell.names,subset.raw=T)
writeMM(mdat@raw.data,file="./Merged1/matrix.mtx")
write.table(cbind(colnames(mdat@raw.data),as.character(mdat@meta.data$group)),file="./Merged1/barcodes.tsv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(cbind(rep(111,length(rownames(mdat@raw.data))),rownames(mdat@raw.data)),file="./Merged1/genes.tsv",col.names=T,row.names=F,sep="\t",quote=F)
