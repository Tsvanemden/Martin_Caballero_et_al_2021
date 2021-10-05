# Authors: Thomas van Emden and Lucia Martin Caballero
# Date: 26.11.2020
# Purpose: Provide extra funcions to reate all NGS figures for MartinCaballero et al., 2021

# Load necessary packages
library(plyr)
library(dplyr)
library(tximport)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(ggrepel)
library(reshape2)
library(sva)
library(BSgenome)
theme_set(theme_pubr())
theme_update(text=element_text(family="Arial"))

# Wrapper for DESeq2 workflow
DESeqFunction <- function(samples, contrast, mutant){
 
  # Some samples are spread over multiple batches, some aren't so this decides if batch correction in needed
  if (length(unique(samples$batch)) >= 2){
    DEseqDesign <- ~ batch + condition
  } else if (length(unique(samples$batch)) == 1){
    DEseqDesign <- ~ condition
  }
  
  # Import the RSEM output
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
                                     design = DEseqDesign)
  
  # Filter genes with low counts
  nrow(ddsTxi)
  keep <- rowSums(counts(ddsTxi)) > 5 * ncol(samples)
  ddsTxi <- ddsTxi[keep,]
  nrow(ddsTxi)
  
  # Do the DESeq
  dds <- DESeq(ddsTxi)
  
  # Prepare the output
  res.df <- toPlot1a <- as.data.frame(results(dds, contrast = contrast))
  colnames(res.df) <- paste(colnames(res.df), mutant, sep = ".")
  res.df$genes <- rownames(res.df)
  
  return(res.df)
}

# Wrapper to produce volcano plots
volcanoPlotter <- function(df, mutant) {
  
  # Prepare df to plot
  toPlot1aSigUp <- subset(df, df$log2FoldChange >= 1 & df$padj < 0.01)
  toPlot1aSigDown <- subset(df, df$log2FoldChange <= -1 & df$padj < 0.01)
  toPlot1aGrey <- df[!(rownames(df) %in% c(rownames(toPlot1aSigUp), rownames(toPlot1aSigDown))),]
  toLabel <- subset(df, df$padj <= 5.389398e-30 | df$log2FoldChange >= 5 )
  
  # Plot data
  volcano <- ggplot(toPlot1aGrey, aes(x=log2FoldChange, y=-log10(padj))) + 
    geom_point(color="grey75", size = 0.25) +
    geom_point(data = toPlot1aSigUp, aes(x=log2FoldChange, y=-log10(padj)), color = "red", size = 0.25) +
    geom_point(data = toPlot1aSigDown, aes(x=log2FoldChange, y=-log10(padj)), color = "blue", size = 0.25) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray50"  ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = -1, linetype = "dashed", color = "gray50") +
    geom_text_repel(data=toLabel, 
                    aes(x=log2FoldChange, y=-log10(padj)),
                    label = rownames(toLabel),
                    colour = "black", size = 1,
                    segment.size = 0.2) +
    xlab(paste0("log2 fold change (WT vs ", mutant, ")")) +
    ylab("-log10(p.adj") +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5))
  
  return(volcano)
}

# Wrapper function to produce single sample IGV style plots
igvPlotterI <- function(data, chr, start, end, fillColor, ylimit) {
  
  # Adjust line number to chromosomal coordinate
  if(chr == "III"){
    nlines <- end - start + 1
    start <- start + 59559 - 1
  } else if (chr == "II"){
    nlines <- end - start + 1
    start <- start + 2512442 - 1
  } else if ( chr == "I"){
    nlines <- end - start + 1
    start <- start + 7052246 - 1
  }
  
  # Read in the files
  toPlot <- read.csv(data, header = FALSE, sep = "\t", nrows = nlines, skip = start)
  
  # Do the plotting
  igv <- ggplot(toPlot, aes(x=V2, y=V3)) +
    geom_area(fill = fillColor) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,ylimit)) +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.ticks = element_blank(),
          axis.line.x = element_blank(),
          plot.margin = unit(c(0,0,0,0), "lines")
    )
  
  return(igv)
}

# Wrapper function to produce 6 sample IGV style plots
igvPlotterIV <- function(df1, df2, df3, df4, df5, df6, chr, start, end, ylimit){
  
  # combine 6 igvPlots
  fig <- ggarrange(igvPlotterI(df1, chr, start, end, "black", ylimit), 
                   igvPlotterI(df2, chr, start, end, "black", ylimit),
                   igvPlotterI(df3, chr, start, end, "black", ylimit), 
                   igvPlotterI(df4, chr, start, end, "red", ylimit),
                   igvPlotterI(df5, chr, start, end, "red", ylimit), 
                   igvPlotterI(df6, chr, start, end, "red", ylimit),
                   nrow = 6 )
  
  return(fig)
}

# Wrapper function to produce pie chart
circlePlotter <- function(dataframe, title){

  # Prepare data by counting element in the df
  circleData <- data.frame(
  group = c("LTR", "ncRNA", "protein coding", "other"),
  value = c(length(grep( "LTR", rownames(dataframe))), 
            length(grep( "NCRNA", rownames(dataframe))), 
            length(grep( "SPAC|SPBC|SPCC", rownames(dataframe))),
            length(rownames(dataframe)) - length(!grep( "SPAC|SPBC|SPCC|NCRNA|LTR", rownames(dataframe))))
  )

  # Set the order of appearance for each element
  circleData$group <- factor(circleData$group , levels = c("ncRNA","protein coding", "LTR", "other"))

  # Plot the pie chart
  piechart <- ggplot(circleData, aes(x="", y=value, fill=group))+
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0, direction = -1) +
    ggtitle(title) +
    scale_fill_manual(values=c("tomato2", "steelblue2", "gray75", "gray25")) +
    ylab(paste0("Total = ", nrow(dataframe))) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 8),
          plot.title = element_text(size = 8, hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.5,"line")
  )
  
  return(piechart)
}

# Wrapper fucntion to produce scatter plots
scatterPlotter <- function(df1, df2, colour){
  
  # Gather all the data from DESeq2 outputs
  toPlot <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "genes"), list(df1, df2))
  
  # Prepare data
  rownames(toPlot) <- toPlot$genes
  toPlot <- subset(toPlot, select=-c(genes))
  toPlot <- toPlot[, grep("log2FoldChange", names(toPlot), value=T)]
  toPlot <- na.omit(toPlot)
  
  # Plot the scatter
  scatterplot <- ggplot(toPlot, aes_string(colnames(toPlot)[1], colnames(toPlot)[2])) + 
    stat_cor(aes(label = ..r.label..), size = 2) +
    geom_point(size = 0.25, color = colour, alpha = 0.5) +
    ylim(min(toPlot), max(toPlot)) +
    xlim(min(toPlot), max(toPlot)) +
    coord_fixed() +
    geom_smooth(method=lm, se=FALSE, alpha=0.5, formula = y ~ x,
                color="black", size = 0.5, linetype="dashed") +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 6)
          )
  
  return(scatterplot)
}