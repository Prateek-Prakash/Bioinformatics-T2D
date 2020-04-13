# Project :: Prateek Prakash (EZ7543)

# Change Working Directory
setwd("~/Desktop/Project")

# Load Libraries
library(Biobase)
library(GEOquery)
library(limma)

# Load Series & Platform Data
gset <- getGEO("GSE20966", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL1352", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Make Proper Column Names To Match Toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Group Names For All Samples
gsms <- "11111111110000000000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# Log2 Transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# Set Up Data & Proceed With Analysis
sml <- paste("G", sml, sep="")
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf, confint = TRUE)

# Boxplot For Select Samples
library(Biobase)
library(GEOquery)

# Load Series & Platform Data
gset <- getGEO("GSE20966", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL1352", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Group Names For All Samples
gsms <- "11111111110000000000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")

# Order Samples By Group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("T2D","Control")

# Set Parameters & Draw Plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1), bg='white')
title <- paste ("GSE20966", '/', annotation(gset), " Selected Samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

dev.copy(png,filename="Exp-Norm.png")
dev.off()

# Annotate Gene Symbols
library(u133x3p.db)
library(annotate)
Gene.Symbols <- getSYMBOL(row.names(tT), "u133x3p")
results <- cbind(Gene.Symbols, tT)

# Data Results (Order By CI.L)
View(results[order(results$CI.L, decreasing = TRUE),])

# Data Results (Order By P.Value)
View(results[order(results$P.Value, decreasing = FALSE),])

# Hierarchical Clustering (RAW Data)
library(rafalib)
clusters = hclust(dist(t(ex)))
myplclust(clusters, lab.col="black", cex=0.8)

dev.copy(png,filename="HC-RAW.png")
dev.off()

# Hierarchical Clustering (Diff. Exp. Data :: CI.L)
library(rafalib)
probes = rownames(results[which(results$CI.L > 1.2),])
clusters = hclust(dist(t(ex[probes,])))
myplclust(clusters, lab.col="black", cex=0.8)

dev.copy(png,filename="HC-Diff-Exp-CI.png")
dev.off()

# Hierarchical Clustering (Diff. Exp. Data :: P.Value)
library(rafalib)
probes = rownames(results[which(results$P.Value < 0.05),])
clusters = hclust(dist(t(ex[probes,])))
myplclust(clusters, lab.col="black", cex=0.8)

dev.copy(png,filename="HC-Diff-Exp-P.png")
dev.off()

# Volcano Plot (Export Manually)
with(tT, plot(logFC, -log10(P.Value), pch=20, main="Volcano Plot", xlim=c(-2.5,2)))

abline(v = c(-1, 1), untf = FALSE)
abline(h = 2, untf = FALSE)

with(subset(tT, P.Value<.01 & abs(logFC)<1), points(logFC, -log10(P.Value), pch=20, col="grey"))
with(subset(tT, P.Value<.01 & logFC>1), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(tT, P.Value<.01 & logFC<(-1)), points(logFC, -log10(P.Value), pch=20, col="green"))

grey.count = nrow(subset(tT, P.Value<.01 & abs(logFC)<1))
grey.count = paste(c("N =", grey.count), collapse = " ")
red.count = nrow(subset(tT, P.Value<.01 & logFC>1))
red.count = paste(c("N =", red.count), collapse = " ")
green.count = nrow(subset(tT, P.Value<.01 & logFC<(-1)))
green.count = paste(c("N =", green.count), collapse = " ")
left.black.count = nrow(subset(tT, P.Value>=.01 & logFC<(-1)))
left.black.count = paste(c("N =", left.black.count), collapse = " ")
right.black.count = nrow(subset(tT, P.Value>=.01 & logFC>1))
right.black.count = paste(c("N =", right.black.count), collapse = " ")

legend(x = "topleft", legend = green.count, bty = "n", text.col = "red")
legend(x = "top", legend = grey.count, bty = "n", text.col = "red")
legend(x = "topright", legend = red.count, bty = "n", text.col = "red")
legend(x = "bottomleft", legend = left.black.count, bty = "n", text.col = "red")
legend(x = "bottomright", legend = right.black.count, bty = "n", text.col = "red")

# Export List Of Genes Of Interest
red.goi = as.character(subset(results, P.Value<.01 & logFC>1)$Gene.Symbols)
red.goi= red.goi[!is.na(red.goi)]
green.goi = as.character(subset(results, P.Value<.01 & logFC<(-1))$Gene.Symbols)
green.goi= green.goi[!is.na(green.goi)]
all.goi = c(red.goi, green.goi)
all.goi = unique(all.goi)
lapply(all.goi, write, "Genes-Of-Interest.txt", append=TRUE, ncolumns=1000)
