# Authors: Thomas van Emden and Lucia Martin Caballero
# Date: 26.11.2020
# Purpose: Create all NGS based figures for Martin Caballero et al., 2021

# Clean environment
rm(list = ls())


# Load extra functions and packages
source("functions_Martin _Caballero_et_al_2021.R")


# Create folder for figures
dir.create("./figures", showWarnings = FALSE )


#### DESeq2 Magic ####


# lem2
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, number >= 19 & strainBackground == 65) 
samples <- subset(samples, condition == "lem2" | condition == "wt")
samples$condition <- relevel(samples$condition, ref = "wt")
res.lem2 <- DESeqFunction(samples = samples, 
                          contrast =c("condition", "lem2" ,"wt"),
                          mutant = "lem2")


# clr4
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, strainBackground == "374")
samples <- subset(samples, condition == "clr4" | condition == "wt") # clr4 and man1 are done in different background so they have different wt
samples$condition <- relevel(samples$condition, ref = "wt")
res.clr4 <- DESeqFunction(samples = samples, 
                          contrast =c("condition", "clr4" ,"wt"),
                          mutant = "clr4")


# man1
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, strainBackground == "374")
samples <- subset(samples, condition == "man1" | condition == "wt") # clr4 and man1 are done in different background so they have different wt
samples$condition <- relevel(samples$condition, ref = "wt")
res.man1 <- DESeqFunction(samples = samples, 
                          contrast =c("condition", "man1" ,"wt"),
                          mutant = "man1")


# rrp6
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, batch == "c")
samples <- subset(samples, condition == "rrp6" | condition == "wt")
samples$condition <- relevel(samples$condition, ref = "wt")
res.rrp6 <- DESeqFunction(samples = samples, 
                          contrast =c("condition", "rrp6" ,"wt"),
                          mutant = "rrp6")


# red1
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, batch == "a")
samples <- subset(samples, condition == "red1" | condition == "wt")
samples$condition <- relevel(samples$condition, ref = "wt")
res.red1 <- DESeqFunction(samples = samples, 
                          contrast =c("condition", "red1" ,"wt"),
                          mutant = "red1")


# erh1
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, batch == "b")
samples <- subset(samples, condition == "erh1" | condition == "wt")
samples$condition <- relevel(samples$condition, ref = "wt")
res.erh1 <- DESeqFunction(samples = samples, 
                          contrast =c("condition", "erh1" ,"wt"),
                          mutant = "erh1")


# air1
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, batch == "e")
samples <- subset(samples, condition == "air1" | condition == "wt")
samples$condition <- relevel(samples$condition, ref = "wt")
res.air1 <- DESeqFunction(samples = samples, 
                          contrast =c("condition", "air1" ,"wt"),
                          mutant = "air1")


# ccr4
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, batch == "b")
samples <- subset(samples, condition == "ccr4" | condition == "wt")
samples$condition <- relevel(samples$condition, ref = "wt")
res.ccr4 <- DESeqFunction(samples = samples, 
                          contrast =c("condition", "ccr4" ,"wt"),
                          mutant = "ccr4")


# iss10
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, batch == "h")
samples <- subset(samples, condition == "iss10" | condition == "wt")
samples$condition <- relevel(samples$condition, ref = "wt")
res.iss10 <- DESeqFunction(samples = samples, 
                          contrast =c("condition", "iss10" ,"wt"),
                          mutant = "iss10")


# lem2iss10
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, batch == "h")
samples <- subset(samples, condition == "lem2iss10" | condition == "wt")
samples$condition <- relevel(samples$condition, ref = "wt")
res.lem2iss10 <- DESeqFunction(samples = samples, 
                           contrast =c("condition", "lem2iss10" ,"wt"),
                           mutant = "lem2iss10")


# Prepare a file that includes all results and export it
res.all <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "genes"), list(res.lem2, res.rrp6, res.red1, res.erh1, res.air1, res.ccr4, res.iss10, res.clr4, res.man1, res.lem2iss10))
rownames(res.all) <- res.all$genes
res.all <- subset(res.all, select=-c(genes))
write.csv(res.all, "figures/resAll.csv")



#### Figure 1 ####
#'#####################################################################
#'##       ###      #####     ##########  #############################
#'##  ##########  ######  ############    #############################
#'##     #######  ######  ##############  #############################
#'##  ##########  ######  ##   #########  #############################
#'##  ##########  ######  ###  #########  #############################
#'##  ########      #####    ########        ##########################
#'#####################################################################


# Fig 1a: Volcano plot of lem2 RNA-seq

# Import and plot the data
fig1Data <- res.lem2
colnames(fig1Data) <- sub(".lem2","",colnames(fig1Data))
fig1a <- volcanoPlotter(df = fig1Data, mutant = "lem2" )


# Fig 1b : MA plot for lem2∆
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, number >= 19 & strainBackground == 65 ) 
samples <- subset(samples, condition == "lem2" | condition == "wt")
samples$condition <- relevel(samples$condition, ref = "wt")

files <- file.path("data/RSEM", paste0(samples$sampleName, ".genes.results"))
names(files) <- samples$sampleName
txi.rsem <- tximport(files, type = "rsem", txIn=FALSE)

# prepare colData
colData <- select(samples, "batch", "condition")
row.names(colData) <- samples$sampleName

# one gene somehow has a length of 0 bp, we have to change this to 1 for the analysis to work
txi.rsem$length[txi.rsem$length == 0] <- 1

# Make the DESeq data set
ddsTxi <- DESeqDataSetFromTximport(txi.rsem,
                                   colData = colData,
                                   design = ~ condition)

# Filter genes with low counts
nrow(ddsTxi)
keep <- rowSums(counts(ddsTxi)) > 5 * ncol(samples)
ddsTxi <- ddsTxi[keep,]
nrow(ddsTxi)

# Do the DESeq and generate results
dds <- DESeq(ddsTxi)
toPlot <- results(dds, alpha=0.05, contrast=c("condition","lem2","wt"))

fig1b <- ggmaplot(toPlot,
         fdr = 0.01, fc = 2, size = 0.005,
         palette = c("#E93325", "#488BCA", "#909090"),
         genenames = as.vector(rownames(toPlot)),
         top = 0,
         alpha = 0.6,
         label.select = c("SPNCRNA.103", "SPSNORNA.34"),
         ggtheme = ggpubr::theme_pubr(legend = "none")) +
         theme(axis.title = element_text(size = 8),
               axis.text = element_text(size = 6),
               axis.line = element_line(size = 0.5))


# Fig 1c: IGV plots of selected genes

# Call the igvPlot function for wt and lem2 for sme2
fig1c1 <- igvPlotterIV("data/bedgraph/wt_7.txt",
                       "data/bedgraph/wt_8.txt",
                       "data/bedgraph/wt_9.txt", 
                       "data/bedgraph/lem2_3.txt", 
                       "data/bedgraph/lem2_4.txt", 
                       "data/bedgraph/lem2_5.txt",
                       chr = "II", start = 339346, end = 341753, ylimit = 55
)

# Add annotation the plots
fig1c1 <- annotate_figure(fig1c1,
                         #left = text_grob("meiRNA", color = "black", rot = 90, size = 8),
                         top = text_grob("sme2+", color = "black", size = 8),
                         fig.lab = "[II:339346-341753]",
                         fig.lab.pos = "top",
                         fig.lab.size = 6) #+ theme(plot.margin = unit(c(0,1.5,0,1.5), "lines"))

# Call the igvPlot function for wt and lem2 for sno20
fig1c2 <- igvPlotterIV("data/bedgraph/wt_7.txt",
                       "data/bedgraph/wt_8.txt",
                       "data/bedgraph/wt_9.txt",
                       "data/bedgraph/lem2_3.txt", 
                       "data/bedgraph/lem2_4.txt", 
                       "data/bedgraph/lem2_5.txt", 
                       chr = "II", start = 2410996, end = 2412526, ylimit = 120
)

# Add annotation the plots
fig1c2 <- annotate_figure(fig1c2,
                          #left = text_grob("small nucleolar RNA", color = "black", rot = 90, size = 8),
                          top = text_grob("sno20+", color = "black", size = 8),
                          fig.lab = "[II:2410996-2412526]",
                          fig.lab.pos = "top",
                          fig.lab.size = 6) #+ theme(plot.margin = unit(c(0,1.5,0,1.5), "lines"))

# Assemble fig1c
fig1c <- ggarrange(fig1c1, fig1c2,
                    ncol = 2, nrow = 1)
                         

# Fig 1d Piecharts of WT vs lem2 upregulated genes

# Get data for whole genome
gtf <- read.csv("genomeFiles/sortedAllLTRs.gtf", header = FALSE, sep = "\t")
gtf <- separate(gtf, V9, sep = ';', into = c("a", "b") ) 
gtf <- gtf[!duplicated(gtf[,"b"]),]
gtf$b <- gsub("gene_id ", replacement = "", gtf$b)
gtf$b <- gsub(" ", replacement = "", gtf$b)
rownames(gtf) <- gtf$b
data1d1 <- gtf

# Plot data whole genome
fig1d1 <- circlePlotter(data1d1, "WT genome\ndistribution")

# Get data for lem2∆
data1d2 <- subset(fig1Data, log2FoldChange >= 1 & padj < 0.01 )

# Plot data lem2∆
fig1d2 <- circlePlotter(data1d2, ">2 fold change\nin lem2")

# Put Fig 1d together
fig1d <- ggarrange(fig1d1, fig1d2,
          common.legend = TRUE,
          legend = "right"
          )


# Assemble fig1
fig1 <- ggarrange(fig1a, fig1b, fig1c, fig1d,
                  nrow = 4, ncol = 1,
                  labels = c("a", "b", "c", "d"),
                  heights = c(67,41,68,24))

# Save figure 1
ggexport(filename =  "figures/fig1.pdf", fig1, width = 8.5/2.54, height = 20/2.54, device = "pdf", useDingbats=FALSE)


#### Supplemental Figure 1 ####
#'#####################################################################
#'##       ###      #####     #######     #######  ####################
#'##  ##########  ######  #########   #########    ####################
#'##     #######  ######  ###########   #########  ####################
#'##  ##########  ######  ##   ########   #######  ####################
#'##  ##########  ######  ###  #########   ######  ####################
#'##  ########      #####    #######     #####        #################
#'#####################################################################


# Fig S1a: Volcano plot of clr4 RNA-seq

# Get and plot the data
figS1aData <- res.clr4
colnames(figS1aData) <- sub(".clr4","",colnames(figS1aData))
figS1a <- volcanoPlotter(df = figS1aData, mutant = "clr4" )


# fig S1b: Volcano plot of man1 RNA-seq

# Get and plot the data
figS1bData <- res.man1
colnames(figS1bData) <- sub(".man1","",colnames(figS1bData))
figS1b <- volcanoPlotter(df = figS1bData, mutant = "man1" )


# Assemble figS1ab
figS1ab <- ggarrange(figS1a, figS1b,
                     nrow = 1, ncol = 2,
                     labels = c("A", "B")
                     )


# fig S1c: Pie charts of various mutants

# Plot data whole genome (same as fig1c1)
figS1c1 <- fig1d1

# Get and plot data for rrp6∆
dataS1c2 <- res.rrp6
colnames(dataS1c2) <- sub(".rrp6","",colnames(dataS1c2))
dataS1c2 <- subset(dataS1c2, log2FoldChange >= 1 & padj < 0.01 )
figS1c2 <- circlePlotter(dataS1c2, ">2 fold change\nin rrp6")

# Get and plot data for red1∆
dataS1c3 <- res.red1
colnames(dataS1c3) <- sub(".red1","",colnames(dataS1c3))
dataS1c3 <- subset(dataS1c3, log2FoldChange >= 1 & padj < 0.01 )
figS1c3 <- circlePlotter(dataS1c3, ">2 fold change\nin red1")

# Get and plot data for clr4∆
dataS1c4 <- figS1aData
colnames(dataS1c4) <- sub(".clr4","",colnames(dataS1c4))
dataS1c4 <- subset(dataS1c4, log2FoldChange >= 1 & padj < 0.01 )
figS1c4 <- circlePlotter(dataS1c4, ">2 fold change\nin clr4")

# Get and plot for man1∆
dataS1c5 <- figS1bData
colnames(dataS1c5) <- sub(".man1","",colnames(dataS1c5))
dataS1c5 <- subset(dataS1c5, log2FoldChange >= 1 & padj < 0.01 )
figS1c5 <- circlePlotter(dataS1c5, ">2 fold change\nin man1")

# Assemble figS1c
figS1c <- ggarrange(figS1c1, figS1c2, figS1c3, figS1c4, figS1c5,
                    nrow = 1, ncol = 5,
                    common.legend = TRUE,
                    legend = "right")

# Assemble figS1
figS1 <- ggarrange(figS1ab, figS1c,
                   nrow = 2, ncol = 1,
                   labels = c("", "C")
)


# Save figure S1
ggexport(filename =  "figures/figS1.pdf", figS1, width = 17/2.54, height = 10/2.54, device = "pdf", useDingbats=FALSE)


###### Figure 2 ####
#'#####################################################################
#'##       ###      #####     ########      ###########################
#'##  ##########  ######  ###########  ####  ##########################
#'##     #######  ######  ################  ###########################
#'##  ##########  ######  ##   #########  #############################
#'##  ##########  ######  ###  #######  ###############################
#'##  ########      #####    ########        ##########################
#'#####################################################################

# fig 2b: k-means plot of RNA-seq data

# prepare the data
# fig 2b: k-means plot of RNA-seq data

# prepare the data
res.fig2b <- res.all[, grep("log2FoldChange", names(res.all), value=T)]
res.fig2b <- res.fig2b[,1:7]

# set the settings

nubKs <- 10 # number of Ks

# Do the Kmeans clustering
set.seed(123)
k <- kmeans(na.omit(res.fig2b), centers = nubKs, iter.max = 1000, nstart = 1)

# Prepare for plotting
dfc <- cbind(na.omit(res.fig2b), id=seq(nrow(na.omit(res.fig2b))), cluster=k$cluster)

# Order clusters by mean of means of genes and samples
dfc$average <- rowMeans(dfc[,1:7])
test <- aggregate(dfc[, 10], list(dfc$cluster), mean)
test$Group.1[test$Group.1 == 1:nubKs] <- letters[seq( from = 1, to = nubKs )]
test <- test[order(test$x),]
test$iddsort <- 1:nubKs

# Reorder clusters
dfc$cluster <- mapvalues(dfc$cluster, from = 1:nubKs, to = letters[seq( from = 1, to = nubKs )])
dfc$cluster <- mapvalues(dfc$cluster, from = test$Group.1, to = test$iddsort)
dfc$cluster <- as.numeric(dfc$cluster)
dfc$idsort <- dfc$id[order(dfc$cluster)]
dfc$idsort <- order(dfc$idsort)

# cap values at 3
dfc[,1:7][dfc[,1:7] >= 3] <- 3
dfc[,1:7][dfc[,1:7] <= -3] <- -3

# Fit cluster values on the scale of the data
dfc$cluster <- as.numeric(dfc$cluster)
dfc$cluster <- (dfc$cluster - (nubKs/2) )/ nubKs * 6

# Make data long and determine order of plotting
dfm <- melt(dfc, id.vars=c("id", "idsort"))
dfm$variable <- sub("log2FoldChange.","",dfm$variable)
dfm$variable <- factor(dfm$variable,levels = c("rrp6", "red1", "lem2","erh1" ,"iss10", "air1", "ccr4", "cluster"))

# Fianl prep for plotting
dfm <- dfm[!is.na(dfm$variable),]
dfm$value <- as.numeric(dfm$value)

# Make the figure
fig2b <- ggplot(dfm, aes(x=variable, y=idsort)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient2(
    low = "blue",
    midpoint = 0,
    mid = "white",
    high = "red") +
  ylab("Log2 Fold Change over WT") +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "left",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.key.size = unit( 0.5, "line")
        )

fig2b

figname <- paste0("Kmeans/", "KmeansPlotFig2_" , nubKs , "Means" , ".pdf")
csvname <- paste0("Kmeans/", "KmeansclusterInfo_" , nubKs , "Means" , ".csv")

ggexport(fig2b, filename = figname, width = 15/2.54, height = 20/2.54, device = "pdf", useDingbats=FALSE)
write.csv(dfc, csvname)


# write dataframe that contains info about what cluster each gene belongs to
write.csv(dfc, "figures/clusterInfo.csv")


# fig 2c: scatter plot of lem2 vs rrp6 RNA-seq data

# plot the scatter
fig2c <- scatterPlotter(res.lem2, res.rrp6, "#5C97C6")


# Assemble fig2
fig2 <- ggarrange(fig2b, fig2c,
                  nrow = 2, ncol = 1,
                  heights = c(6,4),
                  widths = c(6,4),
                  labels = c("B", "C")
                  )  

# Save figure 2
ggsave("figures/fig2.pdf", fig2, width = 11/2.54, height = 10/2.54, device = "pdf", useDingbats=FALSE)


#### Supplemental Figure 2 ####
#'#####################################################################
#'##       ###      #####     #######     #####      ##################
#'##  ##########  ######  #########   ########  ####  #################
#'##     #######  ######  ###########   ###########  ##################
#'##  ##########  ######  ##   ########   #######  ####################
#'##  ##########  ######  ###  #########   ####  ######################
#'##  ########      #####    #######     #####        #################
#'#####################################################################


# fig S2a: PCA plot of the RNA-seq libraries

# Get and prepare the data for DESeq2
samples <- read.table("SraRunTable.txt", header = TRUE)
samples <- subset(samples, samples$condition != "lem2iss10")
files <- file.path("data/RSEM", paste0(samples$sampleName, ".genes.results"))
names(files) <- samples$sampleName
txi.rsem <- tximport(files, type = "rsem", txIn=FALSE)
colData <- select(samples, "batch", "condition")
row.names(colData) <- samples$sampleName

# Export a raw counts file from all samples
rawCountsFile <- txi.rsem$counts
write.csv(rawCountsFile, file="figures/counts.csv")

# one gene somehow has a length of 0 bp, we have to change this to 1 for the analysis to work
txi.rsem$length[txi.rsem$length == 0] <- 1

# Make the DESeq data set
dds <- DESeqDataSetFromTximport(txi.rsem,
                                   colData = colData,
                                   design = ~ condition)

# Filter genes with low counts
keep <- rowSums(counts(dds)) > 5 * ncol(samples)
dds <- dds[keep,]

#do whole pipeline for unfiltered reads, without changing the name of the variables
dds <- DESeq(dds)

# Prepare data for PCA plotting
rldf <- vst(dds, blind=FALSE)
batchVar <- colData(dds)$batch
modcombat <- model.matrix(~ condition, data = colData(dds))
BatchCorrect <- ComBat(dat = assay(rldf), batch = batchVar, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
rldfb <- rldf
assay(rldfb) <- BatchCorrect

# Plot PCA plot
pcaData <- plotPCA(rldfb, intgroup=c("condition"), ntop=6000, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
figS2a <- ggplot(pcaData, aes(PC1, PC2, color=condition)) + 
                 geom_point(size=1) +
                 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                 ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
                 coord_fixed() +
                 theme(axis.title = element_text(size = 8),
                       axis.text = element_text(size = 6),
                       legend.title = element_blank(),
                       legend.position = "right",
                       legend.text = element_text(size = 6),
                       legend.key.size = unit(0.5,"line")
                       )


# fig S2b: scatter plot of lem2 vs red1 RNA-seq data 

# plot the scatter
figS2b <- scatterPlotter(res.lem2, res.red1, "#9283BD")

 
# FigS2c: scatter plot of lem2 vs erh1 RNA-seq data 

# plot the scatter
figS2c <- scatterPlotter(res.lem2, res.erh1, "#696969")


# FigS2d: scatter plot of lem2 vs ccr4 RNA-seq data 

# plot the scatter
figS2d <- scatterPlotter(res.lem2, res.ccr4, "#696969")


# FigS2e: scatter plot of lem2 vs air1 RNA-seq data 

# plot the scatter
figS2e <- scatterPlotter(res.lem2, res.air1, "#696969")


# FigS2f: scatter plot of lem2 vs iss10 RNA-seq data 

# plot the scatter
figS2f <- scatterPlotter(res.lem2, res.iss10, "#696969")


# Assemble fig S2b-e
figS2be <- ggarrange(figS2b,figS2c,figS2d, figS2e,
                   nrow = 1, ncol = 4,
                   labels = c("B", "C", "D", "E" )
                   )

# Assemble fig S2
figS2 <- ggarrange(figS2a,figS2be,
                   nrow = 2, ncol = 1,
                   widths = c(1,2),
                   labels = c("A", "")
                   )

figS2

# Save figure S2
ggsave("figures/figS2.pdf", figS2, width = 17/2.54, height = 12/2.54, device = "pdf", useDingbats=FALSE)


###### Figure 6 ####
#'#####################################################################
#'##       ###      #####     ########     ############################
#'##  ##########  ######  ###########  ################################
#'##     #######  ######  ###########  ################################
#'##  ##########  ######  ##   ######       ###########################
#'##  ##########  ######  ###  ######  ###  ###########################
#'##  ########      #####    #########    #############################
#'#####################################################################


# Kmeans plot with lem2, iss10 and the double

# Get the data
res.fig2b <- res.all[, grep("log2FoldChange", names(res.all), value=T)]
res.fig2b <- res.fig2b[,1:10]
res.fig2b <- res.fig2b[ , grepl( "lem2|iss10" , names( res.fig2b ) ) ]

# Do the Kmeans clustering
nubKs <- 8 # number of Ks
set.seed(123)
k <- kmeans(na.omit(res.fig2b), centers = nubKs, iter.max = 1000, nstart = 1)

# Prepare for plotting
dfc <- cbind(na.omit(res.fig2b), id=seq(nrow(na.omit(res.fig2b))), cluster=k$cluster)

# Order clusters by mean of means of genes and samples
dfc$average <- rowMeans(dfc[,1:3])
test <- aggregate(dfc[, 6], list(dfc$cluster), mean)
test$Group.1[test$Group.1 == 1:nubKs] <- letters[seq( from = 1, to = nubKs )]
test <- test[order(test$x),]
test$iddsort <- 1:nubKs

# Reorder clusters
dfc$cluster <- mapvalues(dfc$cluster, from = 1:nubKs, to = letters[seq( from = 1, to = nubKs )])
dfc$cluster <- mapvalues(dfc$cluster, from = test$Group.1, to = test$iddsort)
dfc$cluster <- as.numeric(dfc$cluster)
dfc$idsort <- dfc$id[order(dfc$cluster)]
dfc$idsort <- order(dfc$idsort)

# Export the unclipped data
csvname <- paste0("figures/", "fig6_clusterInfoUnclipped" , ".csv")
write.csv(dfc, csvname)

# save copy for motif search
motifInput <- dfc

# cap values at 3
dfc[,1:3][dfc[,1:3] >= 3] <- 3
dfc[,1:3][dfc[,1:3] <= -3] <- -3

# Fit cluster values on the scale of the data
dfc$cluster <- as.numeric(dfc$cluster)
dfc$cluster <- (dfc$cluster - (nubKs/2) )/ nubKs*6

# Make data long and determine order of plotting
dfm <- melt(dfc, id.vars=c("id", "idsort"))
dfm$variable <- sub("log2FoldChange.","",dfm$variable)
dfm$variable <- factor(dfm$variable,levels = c("iss10", "lem2", "lem2iss10", "cluster"))

# Fianl prep for plotting
dfm <- dfm[!is.na(dfm$variable),]
dfm$value <- as.numeric(dfm$value)

# Make the figure
fig6 <- ggplot(dfm, aes(x=variable, y=idsort)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient2(
    low = "blue",
    midpoint = 0,
    mid = "white",
    high = "red") +
  ylab("Log2 Fold Change over WT") +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "left",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.key.size = unit( 0.5, "line")
  )

# Save figure S2
ggexport(fig6, filename = "figures/fig6pdf", width = 15/2.54, height = 20/2.54, device = "pdf", useDingbats=FALSE)
write.csv(dfc, "figures/fig6_clusterInfo.csv")


# Prepare the .fasta files for motif search

# Get data for whole genome
gtf <- read.csv("genomeFiles/sortedStrandedAllLTRs.gtf", header = FALSE, sep = "\t")
gtf <- subset(gtf, gtf$V3 == "CDS")
gtf <- separate(gtf, V9, sep = ';', into = c("a", "b") ) 
gtf$b <- gsub("gene_id ", replacement = "", gtf$b)
gtf$b <- gsub(" ", replacement = "", gtf$b)

# Select gene names from res.all, a list or the proviously defined motifInput

# From res.all
res.motif <- res.all
res.motif <- subset(res.motif, res.motif$log2FoldChange.lem2 >= 1)
res.motif <- subset(res.motif, res.motif$log2FoldChange.iss10 < 1)
tof <- (rownames(res.motif))

# From a list YOU like
tof <- c("SPAC32A11.01", "SPAC27D7.13c", "SPNCRNA.103", "SPAC57A10.04", "SPCC70.09c", 
         "SPAC13A11.03", "SPBC29A10.14", "SPBC32H8.11", "SPBC1921.04c", "SPBC337.08c", 
         "SPBC2G2.09c", "SPBC1347.12", "SPBC19F8.04c", "SPAP27G11.08c", "SPCC11E10.03", 
         "SPBC646.17c", "SPBP8B7.04", "SPAC23C4.07", "SPAC6C3.05", "SPAC1556.06", 
         "SPCC1393.07c", "SPBC29A10.02", "SPBC216.02", "SPAC458.04c", "SPBC2D10.06", 
         "SPAC17A5.18c", "SPBC582.06c", "SPCC4E9.01c", "SPAC25G10.04c", "SPBC577.05c")

tof <- c("SPBC32H8.11", "SPBC29A10.02", "SPAC27D7.13c", "SPBC29A10.14")

# From the previously defined motifInput
motifInput <- subset(motifInput, motifInput$cluster == 6)
#motifInput <- motifInput[ grep("NCRNA", rownames(motifInput), value = T),]

tof <- (rownames(motifInput))

# Prepare the genes and coordinates
gtf <- subset(gtf, gtf$b %in% tof)
table1 <- aggregate(V4~V1+b+V7,data=gtf,FUN=min)
table2 <- aggregate(V5~b,data=gtf,FUN=max)
gtf <- merge(table1, table2, by= "b")

# Load the genome
SpombeGenome <- readDNAStringSet("genomeFiles/ASM294v2.27.fasta", format="fasta")
names(SpombeGenome) <- c("MT", "AB325691", "MTR", "III", "II", "I")

# to prevent appending, delete the previous file
nameFile <- "figures/outfileMotifcl8CDS.txt"
if (file.exists(nameFile)) {
  #Delete file if it exists
  file.remove(nameFile)
}

for (i in 1:nrow(gtf)){
  # Get the info needed from the sample table
  transcript <- as.character(gtf$b[i])
  chr <- as.character(gtf$V1[i])
  start <- as.numeric(gtf$V4[i])
  end <- as.numeric(gtf$V5[i])
  strand <- as.character(gtf$V7[i])
  range <- GRanges(chr, IRanges(start=start,end=end), strand=strand)
  
  seq <- getSeq(SpombeGenome, range)
  fasta <- paste0(">", transcript, "\n", as.character(seq))
  
  cat(fasta, file = nameFile, append = TRUE, sep = "\n")
}


