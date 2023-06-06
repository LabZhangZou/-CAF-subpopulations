## ----------------    seurat integration   ---------------- 
library(dplyr)
library(Seurat)
library(harmony)

outdir <- '.'
sample <- c('F86','F60','F53','hNP1','hNP5','hNP4','PAS','PA1','PB1')
group <- c('F','F','F','hNP','hNP','hNP','P','P','P')
exp_list <- c('/path/to/rawdata/F86/filtered_feature_bc_matrix','/path/to/rawdata/FM60/filtered_feature_bc_matrix','/path/to/rawdata/F53/filtered_feature_bc_matrix','/path/to/rawdata/hNP1/filtered_feature_bc_matrix','/path/to/rawdata/xiaodedong0812/filtered_feature_bc_matrix','/path/to/rawdata/huangxiaochun0812/filtered_feature_bc_matrix','/path/to/rawdata/A2_S01/filtered_feature_bc_matrix','/path/to/rawdata/A1020_1/filtered_feature_bc_matrix','/path/to/rawdata/B1020_1/filtered_feature_bc_matrix')

for(i in 1:length(exp_list)){
	if(i == 1){
		input <- Read10X(data.dir = exp_list[i])
		data <- CreateSeuratObject(counts = input, project = sample[i], min.cells = 5)
		data <- RenameCells(data,add.cell.id = sample[i])
		data@meta.data$sample <- sample[i]
		data@meta.data$group <- group[i]
		data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
		ngenemax <- unname(quantile(data$nFeature_RNA,0.98))
		umimax <- unname(quantile(data$nCount_RNA,0.98))
		data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax & percent.mt < 50,nCount_RNA > 0 &  nCount_RNA  < umimax)
		data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
	}else{
		input <- Read10X(data.dir = exp_list[i])
		input <- Read10X(data.dir = exp_list[i])
        tmp <- CreateSeuratObject(counts = input, project = sample[i], min.cells = 5)
        tmp <- RenameCells(tmp,add.cell.id = sample[i])
        tmp@meta.data$sample <- sample[i]
        tmp@meta.data$group <- group[i]
        tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
        ngenemax <- unname(quantile(tmp$nFeature_RNA,0.98))
        umimax <- unname(quantile(tmp$nCount_RNA,0.98))
        tmp <- subset(tmp, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax & percent.mt < 50,nCount_RNA > 0 &  nCount_RNA  < umimax)
        tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 10000)
		data <- merge(data,tmp)
	}
	data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(data)
	data <- ScaleData(data, features = all.genes)
	data <- RunPCA(data, features = VariableFeatures(object = data))
	DimPlot(data, reduction = "pca")
	data <- JackStraw(data, num.replicate = 100)
	data <- ScoreJackStraw(data, dims = 1:20)
	pdf(paste(outdir,'/',prefix,'.PCJackStrawPlot.pdf',sep=''))
	JackStrawPlot(data, dims = 1:20)
	dev.off()
	plot2<-ElbowPlot(data)
	CombinePlots(plots = list(plot1, plot2))
	data <- RunHarmony(PRO,"sample", plot_convergence = FALSE)
	data <- FindNeighbors(data, reduction = "harmony",dims = 1:20)
	data <- FindClusters(data, resolution = 1.2)
	data <- RunUMAP(data, reduction = "harmony",dims = 1:20)
	data <- RunTSNE(data, reduction = "harmony",dims = 1:20)
	pdf(paste(outdir,'/',prefix,'.umap.cluster.pdf',sep=''))
	DimPlot(data, reduction = "umap",label = F)
	dev.off()
	pdf(paste(outdir,'/',prefix,'.tsne.cluster.pdf',sep=''))
	DimPlot(data, reduction = "tsne",label = F)
	dev.off()
	
	for (l in levels (data)) {
	    cluster1.markers <- FindMarkers(data, ident.1 = l, min.pct = 0.01,logfc.threshold = 0.01)
	    cluster1.markers <- data.frame (gene = rownames (cluster1.markers))
	    write.table (cluster1.markers,file.path(outdir,paste0(prefix,".cluster.",l,".diffgenes.xls")),row.names = F,col.names = T, sep='\t',quote = F)
	}
}

## annotation
new.cluster.ids <- celltypes
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
pdf(paste(outdir,'/',prefix,'.umap.label.pdf',sep=''))
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()



## ----------------    Clusterprofiler  gokegg ---------------- 
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

genes <- read.table (diffgene,header = T)
data = bitr(genes$gene
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb="org.Hs.eg.db")

ggo <- groupGO(gene = data$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)
ego_ALL <- enrichGO(gene = data$ENTREZID, 
                universe = names(geneList),
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1)

dotplot(ego_MF,title="EnrichmentGO_MF_dot")
barplot(ego_MF, showCategory=20,title="EnrichmentGO_MF")
kk <- enrichKEGG(gene = data$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 1)
dotplot(kk,title="Enrichment KEGG_dot")


## ## ----------------    Monocle2   ---------------- 
library (monocle)
library (Seurat)

rds <- readRDS (inputdata)
data <- as(as.matrix(rds@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = rds@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)    

mcds=cds
disp_table <- dispersionTable(mcds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
mcds <- setOrderingFilter(mcds, unsup_clustering_genes$gene_id)
plot_ordering_genes(mcds)

plot_pc_variance_explained(mcds, return_all = F)
mcds <- reduceDimension(mcds, max_components = 2, num_dim = 6,
                reduction_method = 'tSNE', verbose = T)
mcds <- clusterCells(mcds, num_clusters = 2)
plot_cell_clusters(mcds, 1, 2, color = "CellType",
    markers = markers)

deg.cluster <- FindAllMarkers(scRNA.Osteoclastic)
diff.genes <- subset(deg.cluster,p_val_adj<0.05)$gene
mcds <- setOrderingFilter(mcds, diff.genes)
plot_ordering_genes(mcds)
var.seurat <- VariableFeatures(scRNA.Osteoclastic)
mcds <- setOrderingFilter(mcds, var.genes)
plot_ordering_genes(mcds)
disp_table <- dispersionTable(mcds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mcds <- setOrderingFilter(mcds, disp.genes)
plot_ordering_genes(mcds)
mcds <- reduceDimension(mcds, max_components = 2,
    method = 'DDRTree')
mcds <- orderCells(mcds)
plot_cell_trajectory(mcds, color_by = "seurat_clusters")
plot_cell_trajectory(mcds, color_by = "State")
plot_cell_trajectory(mcds, color_by = "Pseudotime")
plot_cell_trajectory(mcds, color_by = "State") +
  facet_wrap(~State, nrow = 1)
mcds_expressed_genes <-  row.names(subset(fData(mcds),
                                          num_cells_expressed >= 10))
mcds_filtered <- mcds[mcds_expressed_genes,]
my_genes <- row.names(subset(fData(mcds_filtered),
                             gene_short_name %in% c("YWHAB", "GAPDH", "TNNC1")))
cds_subset <- mcds_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")
plot_genes_in_pseudotime(cds_subset, color_by =  "State")
genes <- c("TNNT2", "TNNC1", "CDK1")
p1 <- plot_genes_jitter(mcds[genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mcds[genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mcds[genes,], color_by = "State")

diff_test_res <- differentialGeneTest(mcds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(mcds[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T)
                        


## ## ----------------    infercnv   ---------------- 
### 
library(Seurat)
library(tidyverse)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(ggpubr)
library(AnnoProbe)

rds <- readRDS (inputdata)

dat <- GetAssayData(rds,assay = "RNA",slot = "counts")
dat <- as.data.frame(dat)
geneFile=read.table (geneFile,header = T)
groupFiles=read.table(groupFiles,header = T)


geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]      
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

dat=dat[match(geneInfor[,1], rownames(dat)),] 
rownames(geneInfor) <- geneInfor$SYMBOL   
geneInfor <- geneInfor[,-1]


meta <- subset(rds@meta.data,select = c("celltype"))


infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=geneInfor,
                                    ref_group_names=refcelltypes)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir=outdir, 
                             cluster_by_groups=FALSE,
                             k_obs_group = 8,
                             denoise=T,
                             HMM=T)

infercnv::plot_cnv(infercnv_obj, 
                   plot_chr_scale = T,
                   output_filename = "better_plot",output_format = "pdf")



## ## ----------------    cellcycle   ---------------- 
#!/usr/bin/env Rscript
suppressMessages({
library(ggplot2)
library(reshape2)
library(argparser)
library(Seurat)
library(dplyr)
library(biomaRt)
library(Cairo)
})

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="the rds file")
argv <- add_argument(argv,"--species", help="the species,human or mouse",default="human")
argv <- add_argument(argv,"--prefix", help="the file prefix")
argv <- add_argument(argv,"--subcluster", help="subcluster list ,split by ,",default='all')
argv <- add_argument(argv,"--threshold", help="genes that are expressed in less than 5%(threshold) of total cells were removed",default=0.05)
argv <- add_argument(argv,"--do_regress", help="do regress out cycling effect",default=FALSE)
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- parse_args(argv)

rds <- argv$rds
species <- argv$species
prefix <- argv$prefix
subcluster <- argv$subcluster
threshold <- as.numeric(argv$threshold)
do_regress <- as.logical(argv$do_regress)
outdir <- argv$outdir
dir.create(outdir)

clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")

PRO <- readRDS(rds)
if (subcluster == 'all') { subcluster = levels(PRO)  }
if (grepl(',',subcluster)){
        substr<-unlist(strsplit(subcluster,split=","))
}else{
        substr<-c(subcluster)
}
print(substr)
PRO <- subset(PRO, idents=substr)


# translate human gene name to mouse ortholog gene name
convertHumanGeneList <- function(x){
	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
	genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
	humanx <- unique(genesV2[, 2])
	return(humanx)
}

if (species == "human"){
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes
} else if (species == "mouse"){
        #s.genes <- convertHumanGeneList(cc.genes$s.genes)
	#g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)
	s.genes <- read.table("m.s.genes.txt",stringsAsFactors=F)
	s.genes <- s.genes$V1
	g2m.genes <- read.table("/m.g2m.genes.txt",stringsAsFactors=F)
	g2m.genes <- g2m.genes$V1
}

# fileter genes that are expressed in less than threshold(default:5%) of total cells
threshold <- threshold * dim(PRO@meta.data)[1]
assay_count <- PRO@assays$RNA@counts
g_use <- intersect(c(s.genes,g2m.genes),rownames(assay_count))
assay_count_use <- assay_count[g_use, ]

count <- as.matrix(assay_count_use)
cell_num <- apply(count,1,function(x) sum(x>0) ) #计算每个基因在多少细胞中有表达量
cell_num <- as.data.frame(cell_num)
cell_num$gene_name <- rownames(cell_num)
pass_gene <- cell_num %>% filter(cell_num > threshold)
pass_gene <- pass_gene$gene_name
s.use <- intersect(s.genes,pass_gene)
g2m.use <- intersect(g2m.genes,pass_gene)
gene.use <- rbind(s.use,g2m.use)
write.table(gene.use,file=paste0(outdir,"/",prefix,"_gene_use.txt"),sep="\t",quote=F,col.names=F)

# calculate cycling score
PRO@meta.data$celltype <- as.character(PRO@active.ident)
PRO <- CellCycleScoring(PRO, s.features = s.use, g2m.features = g2m.use, set.ident = TRUE)
meta <- PRO@meta.data
meta$Phase <- as.character(meta$Phase)
meta$Phase[meta$Phase == "G1"] <- "G0/G1"

# barplot
data=table(meta$Phase,meta$old.ident)
write.table(data,file=paste0(outdir,"/",prefix,"_cell_cycle_num.csv"),col.names=NA,sep=",")
data <- as.data.frame(data)
colnames(data) <- c("phase","celltype","prop")
prop <- prop.table(x=table(meta$Phase,meta$old.ident),margin=2)
write.table(prop,file=paste0(outdir,"/",prefix,"_cell_cycle_prop.csv"),col.names=NA,sep=",")

p <- ggplot(data, aes(x=celltype, y=prop)) + theme_classic() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
geom_bar(stat="identity", position="fill", aes(fill=phase)) +
scale_y_continuous(labels = scales::percent)
ggsave(p,file=paste0(outdir,"/",prefix,"_cellcycle.png"))
ggsave(p,file=paste0(outdir,"/",prefix,"_cellcycle.pdf"))

# before regress out cycling effect
if ("umap" %in% names(PRO@reductions)) {
    p1 <- DimPlot(PRO, reduction = "umap")
    p2 <- DimPlot(PRO, reduction = "umap",group.by = "celltype",cols =clustcol)
} else {
    all.genes <- rownames(PRO)
    PRO_before <- ScaleData(PRO, features = all.genes)
    PRO_before <- RunPCA(PRO_before, features = c(s.use, g2m.use))
    PRO_before <- RunUMAP(PRO_before, reduction = "pca", dims = 1:20)
    p1 <- DimPlot(PRO_before, reduction = "umap")
    p2 <- DimPlot(PRO, reduction = "umap",group.by = "celltype",cols =clustcol)
}

CairoPDF(paste(outdir,'/',prefix,'.umap.CellCycle_before.pdf',sep=''))
print(p1)
dev.off()
CairoPNG(paste(outdir,'/',prefix,'.umap.CellCycle_before.png',sep=''))
print(p1)
dev.off()

CairoPDF(paste(outdir,'/',prefix,'.umap.celltype_before.pdf',sep=''))
print(p2)
dev.off()
CairoPNG(paste(outdir,'/',prefix,'.umap.celltype_before.png',sep=''))
print(p2)
dev.off()

# after regress out cycling effect
if (do_regress) {
    all.genes <- rownames(PRO)
    PRO_after <- ScaleData(PRO, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)
    PRO_after <- RunPCA(PRO_after, features =c(s.use,g2m.use))
    PRO_after <- RunUMAP(PRO_after, reduction = "pca", dims = 1:20)

    p2 <- DimPlot(PRO_after,reduction = "umap")
    CairoPDF(paste(outdir,'/',prefix,'.umap.CellCycle_after.pdf',sep=''))
    print(p2)
    dev.off()
    CairoPNG(paste(outdir,'/',prefix,'.umap.CellCycle_after.png',sep=''))
    print(p2)
    dev.off()
}
# saveRDS(PRO,paste(outdir,'/',prefix,'.cell_cycle.rds',sep=''))



## ----------------    dotplot   ---------------- 
#!/usr/bin/env Rscript
library(argparser)
argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="input rds file, usually generated by rename scirpt")
argv <- add_argument(argv,"--gene", help="gene list file contains genes to draw feature plot, left colume contains genes_raw, right colume contain gene_names, separated with tab")
argv <- add_argument(argv,"--celltypes", help="idents use", default = "all")
argv <- add_argument(argv,"--samples", help="samples use", default = "all")
argv <- add_argument(argv,"--prefix", help="prefix, project id")
argv <- add_argument(argv,"--out", help="the output dir")
argv <- parse_args(argv)

suppressMessages({
    library(argparser)
    library(Seurat)
    library(Cairo)
    library(ggplot2)
})
options(stringsAsFactors=FALSE)

##
rdsfile <- argv$rds
outdir <- argv$out
gene_list <- argv$gene
prefix <- argv$prefix
celltypes <- argv$celltypes
samples <- argv$samples

RDS <- readRDS(rdsfile)

if (! dir.exists(outdir)) dir.create(outdir)

if (grepl("txt",gene_list)) {
    gene_table <- read.table(gene_list,header=F)
    marker <- unique(gene_table$V1)
} else {
    marker <- unlist(strsplit(gene_list,split=","))
}
rnames <- c(rownames(RDS),colnames(RDS[[]]))
marker.tmp <- sapply(marker,function(x) {
    rnames[which(tolower(rnames) == tolower(x))]
})
marker <- unlist(marker.tmp)

if (celltypes != "all") {
    idents <- as.character(unlist(strsplit(celltypes, ",")))
    RDS <- subset(RDS, idents = idents)
}
if (samples != "all") {
    sps <- as.character(unlist(strsplit(samples, ",")))
    RDS <- subset(RDS, subset = sample %in% sps)
}
dot.colors = c("grey","#9e0142")
x.intercept <- 0
ct = anno.celltypes = levels(RDS)

p <- DotPlot(object = RDS, features = marker, col.min = -3, col.max = 3, cols = dot.colors) +
geom_vline(xintercept = x.intercept, colour="#000000B2") +coord_flip() +
scale_y_discrete(breaks=ct, labels=anno.celltypes) +
labs(title = "Gene Expression",x="",y="") +
theme(text = element_text(face = "bold", size = 12),
axis.text.x = element_text(hjust = 1, vjust = 1, angle=45),
plot.title = element_text(size = 16, hjust = 0.5),
panel.background = element_rect(colour = "black",size = 1.2, fill = "transparent",linetype=1))

nMarker <- length(marker)
nCelltype <- length(levels(RDS))

pdf_width <- 4 + 0.2 * nCelltype
pdf_height <- 5 + 0.13 * nMarker
png_width <- 400 + 20 * nCelltype
png_height <- 500 + 13 * nMarker

pdf(paste0(outdir,'/',prefix,'.markergene.DotPlot.pdf'), width = pdf_width, height = pdf_height)
print(p)
dev.off()
CairoPNG(paste0(outdir,'/',prefix,'.markergene.DotPlot.png'), width = png_width, height = png_height)
print(p)
dev.off()

### cord flip
p1 <- DotPlot(object = RDS, features = marker, col.min = -3, col.max = 3, cols = dot.colors) +
geom_vline(xintercept = x.intercept, colour="#000000B2") +
scale_y_discrete(breaks=ct, labels=anno.celltypes) +
labs(title = "Gene Expression",x="",y="") +
theme(text = element_text(face = "bold", size = 12),
axis.text.x = element_text(hjust = 1, vjust = 1, angle=45),
plot.title = element_text(size = 16, hjust = 0.5),
panel.background = element_rect(colour = "black",size = 1.2, fill = "transparent",linetype=1))

nMarker <- length(marker)
nCelltype <- length(levels(RDS))

pdf_height <- 3 + 0.1 * nCelltype
pdf_width <- 5 + 0.2 * nMarker
png_height <- 200 + 10 * nCelltype
png_width <- 500 + 13 * nMarker

pdf(paste0(outdir,'/',prefix,'.markergene.DotPlot.flip.pdf'), width = pdf_width, height = pdf_height)
print(p1)
dev.off()
CairoPNG(paste0(outdir,'/',prefix,'.markergene.DotPlot.flip.png'), width = png_width, height = png_height)
print(p1)
dev.off()


## ----------------    feature plot   ---------------- 
#!/usr/bin/env Rscript
library(argparser)
argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="input rds file, usually generated by rename scirpt")
argv <- add_argument(argv,"--gene", help="gene list file contains genes to draw feature plot, left colume contains genes_raw, right colume contain gene_names, separated with tab")
argv <- add_argument(argv,"--prefix", help="prefix, project id")
argv <- add_argument(argv,"--out", help="the output dir")
argv <- parse_args(argv)

suppressMessages({
	library(argparser)
	library(Seurat)
	library(Cairo)
})
options(stringsAsFactors=FALSE)

##
rdsfile <- argv$rds
outdir <- argv$out
gene_list <- argv$gene
prefix <- argv$prefix

RDS <- readRDS(rdsfile)

if (! dir.exists(outdir)) dir.create(outdir)
if (grepl("txt",gene_list)) {
    gene_table <- read.table(gene_list,header=F)
    marker <- unique(gene_table$V1)
} else {
    marker <- unlist(strsplit(gene_list,split=","))
}

rnames <- c(rownames(RDS),colnames(RDS[[]]))
marker.tmp <- sapply(marker,function(x) {
    rnames[which(tolower(rnames) == tolower(x))]
})

marker <- unlist(marker.tmp)

clust<-summary(RDS@active.ident)
cluster_cell<-as.data.frame(clust)
cell_number<-sum(cluster_cell$clust)
print(cell_number)
pt_use<-0.6
if(cell_number>1000){
	pt_use<-0.4
}
if(cell_number>2500){
	pt_use<-0.3
}
if(cell_number>4000){
	pt_use<-0.2
}
if(cell_number>5500){
	pt_use<-0.15
}
if(cell_number>6500){
	pt_use<-0.1
}


glen <- length(marker)
pload <- 4
pdf_width = 7
pdf_height = 7
png_width = 480
png_height = 480

if (pload == 4){
	pdf_width = 8
	pdf_height = 7
    png_width = 580
    png_height = 480
}

div <- glen %/% pload
mod <-  glen %% pload
pnum <- 0
if (div > 0) {
	for (i in 1:div) {
		g_start <- pload * (i-1) + 1
		g_end <- pload * i
		plot_g<- marker[g_start:g_end]
        pnum.name <- paste0(".",paste(plot_g,collapse="."))
		pnum <- pnum + 1
		p1 <- FeaturePlot(RDS,features=plot_g,cols = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use,order=T,reduction="tsne")
		pdf(paste0(outdir,'/',prefix,'.markergene_tsneFeatureplot',pnum.name,'.pdf'), width = pdf_width, height = pdf_height)
		print(p1)
		dev.off()
		CairoPNG(paste0(outdir,'/',prefix,'.markergene_tsneFeatureplot',pnum.name,'.png'), width = png_width, height = png_height)
		print(p1)
		dev.off()
		p2 <- FeaturePlot(RDS,features=plot_g,cols = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use, order=T,reduction="umap")
		pdf(paste0(outdir,'/',prefix,'.markergene_umapFeatureplot',pnum.name,'.pdf'), width = pdf_width, height = pdf_height)
		print(p2)
		dev.off()
		CairoPNG(paste0(outdir,'/',prefix,'.markergene_umapFeatureplot',pnum.name,'.png'), width = png_width, height = png_height)
		print(p2)
		dev.off()
	}
}

if (mod > 0 ) {
	if (mod == 2){
		pdf_width = 8
		pdf_height = 3.5
        png_width = 580
        png_height = 240
	} else if (mod == 1) {
		pdf_width = 3.5
		pdf_height = 3.5
        png_width = 240
        png_height = 240
	}
	g_start <- pload * div + 1
	g_end <- pload * div + mod
	plot_g<- marker[g_start:g_end]
    pnum.name <- paste0(".",paste(plot_g,collapse="."))
	pnum <- pnum + 1
	p1 <- FeaturePlot(RDS,features=plot_g,cols = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use, order=T,reduction="tsne")
	pdf(paste0(outdir,'/',prefix,'.markergene_tsneFeatureplot',pnum.name,'.pdf') , width = pdf_width, height = pdf_height)
	print(p1)
	dev.off()
	CairoPNG(paste0(outdir,'/',prefix,'.markergene_tsneFeatureplot',pnum.name,'.png'), width = png_width, height = png_height)
	print(p1)
	dev.off()
	p2 <- FeaturePlot(RDS,features=plot_g,cols = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use, order=T,reduction="umap")
	pdf(paste0(outdir,'/',prefix,'.markergene_umapFeatureplot',pnum.name,'.pdf'), width = pdf_width, height = pdf_height)
	print(p2)
	dev.off()
	CairoPNG(paste0(outdir,'/',prefix,'.markergene_umapFeatureplot',pnum.name,'.png'), width = png_width, height = png_height)
	print(p2)
	dev.off()

}




## ----------------    violin plot   ---------------- 
#!/usr/bin/env Rscript
suppressMessages({
	library(argparser)
	library(Seurat)
	library(Cairo)
	library(ggplot2)
})
options(stringsAsFactors=FALSE)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="input rds file, usually generated by rename scirpt")
argv <- add_argument(argv,"--gene_list", help="gene list file contains genes to draw feature plot, left colume contains genes_raw, right colume contain gene_names, separated with tab")
argv <- add_argument(argv,"--prefix", help="prefix, project id")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--VlnDotSize", help="dot size of vln plot", default = 0.01)
argv <- add_argument(argv,"--newcol", help="using new color protocol", default='T')
argv <- parse_args(argv)
##
rdsfile <- argv$rds
outdir <- argv$outdir
gene_list <- argv$gene_list
prefix <- argv$prefix
VlnDotSize <- as.numeric(argv$VlnDotSize)
## stat
clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")

if(argv$newcol=='T'){
source('color_protocol.R')
clustcol <- c(color_protocol, clustcol)
}


RDS <- readRDS(rdsfile)

if (! dir.exists(outdir)) dir.create(outdir)
gene_table <- read.table(gene_list,header=F,stringsAsFactors=F)
if(gene_table$V1[1] == TRUE){
	gene_table$V1[1] = 'T'
}
marker <- unique(gene_table$V1)
marker <- intersect(marker, c(rownames(RDS), colnames(RDS@meta.data)))

glen <- length(marker)
pload <- 4

div <- glen %/% pload
mod <-  glen %% pload
pnum <- 0
if (div > 0) {
	for (i in 1:div) {
		g_start <- pload * (i-1) + 1
		g_end <- pload * i
		plot_g<- marker[g_start:g_end]
		pnum <- pnum + 1
		gname <- paste(plot_g, collapse='_')
		pdf(paste0(outdir,'/',prefix, '.markergene_vlnplot',pnum,'.',gname,'.pdf'), width = 10)
		print(VlnPlot(RDS,features=plot_g,pt.size= VlnDotSize, cols=clustcol,ncol=2))
		dev.off()
		CairoPNG(paste0(outdir,'/',prefix, '.markergene_vlnplot',pnum,'.',gname,'.png'), width = 1000)
		print(VlnPlot(RDS,features=plot_g,pt.size= VlnDotSize, cols=clustcol,ncol=2))
		dev.off()
	}
}

if (mod > 0 ) {
	g_start <- pload * div + 1
	g_end <- pload * div + mod
	plot_g<- marker[g_start:g_end]
	pnum <- pnum + 1
	gname <- paste(plot_g, collapse='_')
	if (mod == 1) {
		pdf_width = 7
		pdf_height = 3.5
		png_width = 700
		png_height= 240
	}
	if (mod == 2) {
		pdf_width = 10
		pdf_height = 3.5
		png_width = 1000
		png_height= 240
	}
	if (mod == 3) {
		pdf_width = 10
		pdf_height = 7
		png_width = 1000
		png_height= 480
	}
	p = VlnPlot(RDS,features=plot_g,pt.size= VlnDotSize, cols=clustcol,ncol=2) +
		theme(axis.text=element_text(size=7),axis.title=element_text(size=10))
	pdf(paste0(outdir,'/',prefix, '.markergene_vlnplot',pnum,'.',gname,'.pdf'), width = pdf_width, height = pdf_height)
	print(p)
	dev.off()
	CairoPNG(paste0(outdir,'/',prefix, '.markergene_vlnplot',pnum,'.',gname,'.png'), width = png_width, height = png_height)
	print(p)
	dev.off()
	p = VlnPlot(RDS,features=plot_g,pt.size= 0.01, cols=clustcol,ncol=2) +
		theme(axis.text=element_text(size=7),axis.title=element_text(size=10))
	pdf(paste0(outdir,'/',prefix, '.markergene_vlnplot',pnum,'.',gname,'.dot.pdf'), width = pdf_width, height = pdf_height)
	print(p)
	dev.off()
	CairoPNG(paste0(outdir,'/',prefix, '.markergene_vlnplot',pnum,'.',gname,'.dot.png'), width = png_width, height = png_height)
	print(p)
	dev.off()
}



## ----------------    scGSVA heatmap plot   ---------------- 
library(dplyr)
library(ggplot2)
library(tidyverse)
library(scGSVA)
library(Seurat)
library(pheatmap)
library(viridis)

Clustering = function(x) {
    d <- dist(x, method = "euclidean")
    fit2 <- hclust(d, method="ward.D")
    id = rownames(x)[fit2$order]
    return(id)
}

SubSetbyPercent = function(pro,pct) {
    cells = colnames(pro)
    cell_use = sample(cells, size= round(length(cells) * pct), replace=F)
    cell_use_order = intersect(cells, cell_use)
    pro_sub = subset(pro, cells=cell_use)
    return(pro_sub)
}

outdir = '.'

clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4",
"#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC67
3","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9")

source('color_protocol.R')
clustcol <- c(color_protocol, clustcol)


PRO = readRDS(rdspaht)

score = readRDS('sample_gsva_score.rds')
gsva = score@gsva
meta = PRO[[]]
meta_new = cbind(meta, gsva[rownames(gsva),])
PRO@meta.data = meta_new


t_gsva = t(gsva)
id_order = Clustering(t_gsva)

PRO_sub = SubSetbyPercent(PRO, 0.3)
t_gsva_scaled = scale(t_gsva)

scaleddata = t_gsva_scaled[,colnames(PRO_sub)]


PRO_sub@assays$RNA@scale.data = scaleddata


p <- DoHeatmap(PRO_sub, features = id_order, label = T, raster = F, draw.lines = F,group.colors = clustcol) +
        scale_fill_gradient2(low = "#42618F", mid = "white", high = "#B43232") +
        theme(axis.text.y = element_text(size = rel(0.8), hjust = 1, colour = "black"))
  
png(file.path(outdir, 'P21100701_No.hallmark.gsva.heatmap.png'))
print(p)
dev.off()
pdf(file.path(outdir, 'P21100701_No.hallmark.gsva.heatmap.pdf'))
print(p)
dev.off()


## topn 
heatmap_data = read.delim('mean_heatmap_data.xls', header =T, row.names=1)
id_use = rownames(heatmap_data)


p <- DoHeatmap(PRO_sub, features = id_use, label = T, raster = F, draw.lines = F,group.colors = clustcol) +
        scale_fill_gradient2(low = "#42618F", mid = "white", high = "#B43232") +
        theme(axis.text.y = element_text(size = rel(0.8), hjust = 1, colour = "black"))
  
png(file.path(outdir, 'P21100701_No.hallmark.gsva.heatmap.top.png'))
print(p)
dev.off()
pdf(file.path(outdir, 'P21100701_No.hallmark.gsva.heatmap.top.pdf'))
print(p)
dev.off()

### mean heatmap
d3 = data.frame()
celltypes = levels(PRO)
for (n in c(1:length(celltypes))) {
    i = celltypes[n]
    barcode = Idents(PRO)[Idents(PRO) == i] %>% names()
    t_gsva_sub = apply(t_gsva[,barcode], 1, mean)
    df_sub = data.frame(i=t_gsva_sub)
    if (n == 1) {
        d3 = df_sub
    } else {
        d3 = cbind(d3, df_sub)
    }
}
colnames(d3) = celltypes

p <- pheatmap(d3,
            cluster_rows = T,
            cluster_cols = F,
            fontsize = 8,
            scale = "row",
            color = colorRampPalette(c("#49679B","white","#E7505C"))(100))

png(file.path(outdir, 'P21100701_No.hallmark.gsva.heatmap.mean.png'),width=500)
print(p)
dev.off()
pdf(file.path(outdir, 'P21100701_No.hallmark.gsva.heatmap.mean.pdf'),width=6)
print(p)
dev.off()