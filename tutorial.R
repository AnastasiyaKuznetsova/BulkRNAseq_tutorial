###########################################################################################
#########################          DATA LOADING          ##################################
###########################################################################################

# Install if you don’t have it
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

# Load package
library(GEOquery)
library(edgeR)
library(Homo.sapiens)
library(Glimma)

eset <- getGEO("GSE104288", GSEMatrix = T)[[1]]
expr <- exprs(eset)
head(expr)

phenoData <- pData(eset)
# featureData is an AnnotatedDataFrame that stores metadata about the features (genes) of the experiment.
# exmaples of what could be found in featureData: entrezID, genomic location
# We lack feature Data
featureData <- fData(eset)

# Here is the expression matrix
expr <- read.csv('GSE104288_DESeq_counts.txt', sep='\t') 
head(expr)
# Let's name rownames as gene names
rownames(expr) <- expr$gene
# Let's get rid of the first gene column which is redundant
expr <- expr[-1]

# write GSM identifiers to a separate column
phenoData['GSM'] <- rownames(phenoData)
# make rownames of phenodata sample identifiers
rownames(phenoData) <- phenoData$title
# match phenoData rows and expr colnames
phenoData <- phenoData[colnames(expr),]

# we'll create featureData object by adding a column corresponding to expr rownames
featureData <- data.frame(geneID = rownames(expr))
rownames(featureData) <- rownames(expr)

# Create eset object for Limma - not sure I need it here
eset <- ExpressionSet(assayData = as.matrix(expr),
                      phenoData = AnnotatedDataFrame(phenoData),
                      featureData = AnnotatedDataFrame(featureData))
dim(eset)

df <- read.csv('GSE104288_DESeq_counts.txt', sep='\t') 
head(df) 
group <- factor(phenoData$characteristics_ch1)
# Create DGEList object
dge <- DGEList(counts = df, group = group)

# Inspect
dge
###########################################################################################
######################          DATA PREPROCESSING          ###############################
###########################################################################################

# Let's annotate our genes - feature matrix

# Example: vector of gene IDs (can be Ensembl, Entrez, etc.)
genes <- featureData$geneID

# Retrieve SYMBOL (gene symbol) and CHR (chromosome)
gene_info <- select(Homo.sapiens,
                    keys = as.character(genes),
                    keytype = "SYMBOL",   # change to "ENTREZID" or "SYMBOL" if needed
                    columns = c("ENSEMBL", "GENENAME","CHR"))
head(gene_info,10)
# Keep only first occurences of each gene ENSEMBLE ID
gene_info_unique <- gene_info[!duplicated(gene_info$SYMBOL), ]
rownames(gene_info_unique) <- gene_info_unique$SYMBOL

""" For differential expression and related analyses, gene expression is rarely considered at the level of raw counts 
since libraries sequenced at a greater depth will result in higher counts. Rather, it is common practice to transform 
raw counts onto a scale that accounts for such library size differences.Popular transformations include counts per million (CPM), 
log2-counts per million (log-CPM), reads per kilobase of transcript per million (RPKM), and fragments per kilobase of transcript 
per million (FPKM)."""

"""In our analyses, CPM and log-CPM transformations are used regularly although they do not account for gene length differences as RPKM 
and FPKM values do. Whilst RPKM and FPKM values can just as well be used, CPM and log-CPM values can be calculated using a counts 
matrix alone and will suffice for the type of comparisons we are interested in."""

lcpm <- cpm(dge, log=TRUE)

# Removing genes that are lowly expressed; characteristics_ch1 - disease status
keep.exprs <- filterByExpr(dge, group=phenoData$characteristics_ch1)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)

library(RColorBrewer)
par(mar = c(5, 4, 4, 2) + 0.1)
nsamples <- ncol(dge)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
lcpm <- cpm(dge, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

# Normalising gene expression distributions
dge <- calcNormFactors(dge, method = "TMM")
dge$samples$norm.factors

# Unsupervised clustering of samples
lcpm <- cpm(dge, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="Sample groups")

# Interactive plot
glMDSPlot(lcpm, labels=paste(group, sep="_"), 
          groups=dge$samples, launch=T)

###########################################################################################
##############           Differential expression analysis           #######################
###########################################################################################
age <- phenoData$`sampling age:ch1`
age_list <- c()

for (a in age) {
  if (grepl("d",a)) {
    a <- 0
  }
    else {
      a <- as.numeric(substr(a, 1, 2))
    }
      
  age_list <- c(age_list,a)
}

# rename names in group
library(dplyr)
group_new <- recode(group, "genotype: Control" = "Control", "genotype: Friedreich's ataxia" = "Ataxia")

# Create a design matrix, adding age as covariate
design <- model.matrix(~group_new+age_list)
colnames(design) <- gsub("group_new", "", colnames(design))
design

dev.off()
par(mfrow=c(1,2))
v <- voom(dge, design, plot=TRUE)
vfit <- lmFit(v, design)
contrast.mat <- makeContrasts(
  disease_vs_WT = Ataxia,  # 1 * groupdisease
  levels = design
)
fit2 <- contrasts.fit(vfit, contrast.mat)
efit <- eBayes(fit2)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))
tfit <- treat(vfit, lfc=0.2)
dt <- decideTests(tfit,)
summary(dt)

tfit <- treat(vfit, lfc=1, trend=T)
dt <- decideTests(tfit,)
summary(dt)

# See distributions of logfold changes
dev.off()
hist(fit2$coefficients[,"disease_vs_WT"], breaks=50, main="Log2FC distribution")

###########################################################################################
#########################           REFERENCES           ##################################
###########################################################################################
"""
I used these resouces for the analysis
Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104288&utm_source=chatgpt.com 
Tutorial: https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#reading-in-count-data

Check this tutorial: https://scienceparkstudygroup.github.io/rna-seq-lesson/ 
Check this video with the data: https://www.youtube.com/watch?v=KJDnk6Glbgg 
Check this: https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf 
"""