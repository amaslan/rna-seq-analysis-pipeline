# differential expression analysis
# Annie Maslan
# usage:
# Rscript --vanilla de.R <output directory> <csv file with sample-cond> 
# references:
# https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html
# https://www.bioconductor.org/packages/3.7/bioc/vignettes/limma/inst/doc/usersguide.pdf

library(ensembldb)
library(tximport)
library(edgeR)
library(limma)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Must provide 2 arguments: (i) output directory (ii) csv file of samples-conditions", call.=FALSE)
}

dir = args[1]
sampleKey = read.csv(file=args[2], header=TRUE, sep=",")
samples = sampleKey$sample
conds = sampleKey$cond

files <- file.path(dir, samples, "quant.sf")
names(files) <- samples 

###########################################################
# CONVERT FROM COUNTS BY TRANSCRIPT TO COUNTS BY GENE
# get table of transcript id - gene id link with gene name as well

# wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz 
# --> from: http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/ --> for de script
gtffile <- "Homo_sapiens.GRCh38.104.gtf.gz"
DB <- ensDbFromGtf(gtf=paste0(ensemblhost, gtffile))
EDB <- EnsDb(DB)
txdf <- transcripts(EDB,
                    return.type="DataFrame",
                    columns=c('gene_id', 'gene_name'))
tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id", "gene_name")])

# read in salmon counts and summarize at gene level
txi <- tximport(files, type = "salmon", txOut=TRUE, countsFromAbundance="lengthScaledTPM")
countGene <- summarizeToGene(txi, tx2gene, ignoreTxVersion = TRUE)

gene_counts <- as.data.frame(countGene$counts)
gene_counts$gene_id <- rownames(gene_counts)
gene_counts$gene <- tx2gene[match(gene_counts$gene_id, tx2gene$gene_id),3]
write.csv(as.data.frame(gene_counts), file = paste(dir, "counts_by_gene.csv", sep="/"))

###########################################################
# DETERMINE DIFFERENTIAL EXPRESSION

# rownames of sampleTable align with the colnames of countGene$counts
sampleTable <- data.frame(condition = factor(conds))
rownames(sampleTable) <- colnames(countGene$counts)

# countGene$counts: matrix of read counts with rows for genes and columns for samples
# create DGEList object with edgeR
dge <- DGEList(countGene$counts, remove.zeros = TRUE)
dge <- calcNormFactors(dge)

# design matrix which indicates which RNA samples have been applied to each array
design <- model.matrix(~ 0+factor(conds))
colnames(design) <- substr(colnames(design),14,nchar(colnames(design)))

# Transform count data to log2-counts per million (logCPM), estimate the mean-variance 
# relationship and use this to compute appropriate observational-level weights. 
# The data are then ready for linear modelling.
# also create dispersion plot
png(paste(dir, "dispersion.png", sep="/"))
v <- voom(dge, design, plot=TRUE)
dev.off()

# apply limma pipelines for differential expression

# start by fitting a linear model to data which fully models the systematic part of the data
fit <- lmFit(v, design)
# empirical Bayes method to moderate the standard errors of the estimated log-fold changes
fit <- eBayes(fit)

# contrast matrix which specifies which comparisons you would like to make between the RNA samples
# To make all pair-wise comparisons between the groups the appropriate contrast matrix can be created by:
design.pairs <- function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
  }
contr.matrix <- design.pairs(colnames(fit$design))

# fit the contrast matrix to the previous fitted model
fit2 <- contrasts.fit(fit, contr.matrix)
# eBayes will compute the consesus pooled variance and then uses
# it to compute empirical Bayes pooled variance for each gene
fit2 <- eBayes(fit2)


##############################
# VOLCANO & P-VAL PLOTS 

# volcano plot
# A volcano plot displays log fold changes on the x-axis versus a measure of statistical significance
# on the y-axis. Here the significance measure can be -log(p-value) or the B-statistics, which give the
# posterior log-odds of differential expression.
# The plot is optionally annotated with the names of the most significant genes.

# p-value plot
# plot p-values (with null hypothesis true for every comparison and p-values
# would be uniformly distributed and expect false positives)
for (i in 1:ncol(contr.matrix)) {
  png(paste(dir, paste("volcano_", i, ".png", sep=""), sep="/"))
  volcanoplot(fit2, coef = i, cex=0.25, main=colnames(fit2$p.value)[i])
  dev.off()
  
  png(paste(dir, paste("phist_", i, ".png", sep=""), sep="/"))
  hist(fit2$p.value[,i],main=colnames(fit2$p.value)[i], xlab = 'p-value')
  dev.off()
}

##############################

# here we use a stricter definition of significance (may consider not using treat if fewer 
# differences in expression)
# treat method (McCarthy and Smyth 2009) can be used to calculate p-values from empirical 
# Bayes moderated t-statistics with a minimum log-FC requirement. The number of differentially 
# expressed genes are reduced when testing requires genes to have a log-FC that is 
# significantly greater than 1 (equivalent to a 2-fold difference between cell types on the 
# original scale).
# the alternate hypothesis you are testing is that the absolute difference between the 
# mean of two groups is greater than some quantity (rather than the conventional alternate 
# hypothesis, which is that the absolute difference is greater than zero).

tfit <- treat(fit2, lfc=1)
dt <- decideTests(tfit)
summary <- summary(dt)
write.csv(summary, file = paste(dir, "de_summary.csv", sep="/"))

# create list of dataframes with significant differential expression for each pairwise comparison
# write each dataframe to csv
sig_tables <- list()
for (i in 1:ncol(contr.matrix)) {
  sig_tables[[i]] <- topTreat(tfit, coef=i, sort.by = 'P', p.value=0.05, n=Inf)
  sig_tables[[i]]$gene_id <- rownames(sig_tables[[i]])
  sig_tables[[i]]$gene<-tx2gene[match(sig_tables[[i]]$gene_id, tx2gene$gene_id),3]
  write.csv(sig_tables[[i]], file = paste(dir, paste("diff_pairwise_", colnames(contr.matrix)[i], ".csv", sep=""), sep="/"))
}
