# rna_wrangle
My assignments #6
---
title: "Rna_wrangle.rmd"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(dplyr)
setwd("~/RNASeqExample/")
samples <- read.csv('sample_info.csv',header = TRUE, sep = ",",
quote = "\"", dec = ".", fill = TRUE, row.names = 1)
genes <- read.csv('expression_results.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
d <- density(genes$KITA_02) # returns the density data
plot(d)
plot(density(log2(genes$KITA_02[(genes$KITA_02>0)])))
```

```{r}
library(dplyr)
setwd("~/RNASeqExample/")
samples <- read.csv('sample_info.csv',header = TRUE, sep = ",",
quote = "\"", dec = ".", fill = TRUE, row.names = 1)
genes <- read.csv('expression_results.csv',header = TRUE, sep = ",",
quote = "\"", dec = ".", fill = TRUE, row.names = 1)
d <- density(samples$PF_BASES)
plot(d)
```

```{r}
library(dplyr)
setwd("~/RNASeqExample/")
samples <- read.csv('sample_info.csv',header = TRUE, sep = ",",
                    quote = "\"", dec = ".", fill = TRUE, row.names = 1)
genes <- read.csv('expression_results.csv',header = TRUE, sep = ",",
                  quote = "\"", dec = ".", fill = TRUE, row.names = 1)
d <- density(genes$KITA_01)
plot(d)
plot(density(log2(genes$KITA_01[(genes$KITA_01>0)])))
d <- density(genes$KITA_03)
plot(d)
plot(density(log2(genes$KITA_03[(genes$KITA_03>0)])))
plot(log2(genes$KITA_01[(genes$KITA_01>0 |genes$KITA_03>0 )]),log2(genes$KITA_03[(genes$KITA_01>0 |genes$KITA_03>0 )]))
```

```{r}
library(ggplot2)
library(reshape2)
corr<-cor(genes)
melted_corr <- melt(corr)
head(melted_corr)
p<-ggplot(melted_corr , aes(x = Var1, y = Var2)) + geom_raster(aes(fill = value)) + scale_fill_gradient2(low="green" , high="red", midpoint=0.5) + theme( plot.title = element_blank(),axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank())
plot(p)
```

```{r}
library('dendextend')
dend <- as.dendrogram(clusters)
dend <- rotate(dend, 1:93)
dend <- color_branches(dend, k=4)
par(cex=0.5)
plot(dend)
```

```{r}
library(ggplot2)
library(plotly)
setwd("~/RNASeqExample/")
samples <- read.csv('sample_info.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
genes <- read.csv('expression_results.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
min(genes[genes>0])
8.05e-12
genes.log <-log2(genes+8.05e-12)
genes.log.small <- genes.log[seq(1, nrow(genes.log), 20), ]
pca <- prcomp(genes.log.small,center = TRUE,scale. = TRUE)
plot(pca, type = "l")
std_dev <- pca$sdev
pr_var <- std_dev^2
pr_var[1:10]
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex, xlab = "Principal Component", ylab = "Proportion of Variance Explained", type = "b")
pcadf<-data.frame(pca$rotation)
plot_ly(data = pcadf, x = ~PC1, y = ~PC2, text = rownames(pcadf))
plot_ly(pcadf, x = ~PC2, y = ~PC3, z = ~PC5) %>%
 add_markers() %>%
 layout(scene = list(xaxis = list(title = 'PC2'), yaxis = list(title = 'PC3'), zaxis = list(title = 'PC5')))
```
