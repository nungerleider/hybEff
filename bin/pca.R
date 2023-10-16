library(SummarizedExperiment)
library(readr)
library(DESeq2)
library(ggplot2)
library(vsn)
library(DescTools)



# Taken from DESeq2 package and altered to allow choosing PCs 
plotPCA = function(object, intgroup="condition", ntop=500, returnData=FALSE, pca1=1, pca2=2)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pca1], PC2=pca$x[,pca2], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
        coord_fixed()
    #return(pca)
}


PCAdf = function(object, intgroup="condition", ntop=500, returnData=FALSE,pca1=1,pca2=2)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pca1], PC2=pca$x[,pca2], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
        coord_fixed()
    return(pca)
}


# Burkitt's Data

# Clinical table
pheno <- read_csv('/Users/nate/Projects/EBV_interactome/bl/blood871418-suppl2.csv')

df <- data.frame(pheno)
df <- data.frame(lapply(df, factor))
colnames(df) <- c('cohort', 'barcode', 'variant', 'ebvstatus', 'genome','sex', 'age', 'biopsy', 'translocation', 'isotype','site', 'source')
df$barcode <- lapply(df$barcode, function(y) gsub("-",".",y))

# EBV genes were removed - only human genes considered
counts <-read.csv('/Users/nate/Projects/EBV_interactome/bl_no_ebvgenes.tsv', sep='\t')
counts_matrix = counts[,2:113]
rownames(counts_matrix) <- counts[,1]
counts <- as.matrix(round(counts_matrix,0))

se <- SummarizedExperiment(assays=counts)
m <- match(colnames(se), df$barcode)
se$cohort <-df$cohort[m]
se$barcode <-df$barcode[m]
se$variant <-df$variant[m]
se$ebvstatus <-df$ebvstatus[m]
se$genome <-df$genome[m]
se$sex <-df$sex[m]
se$age <-df$age[m]
se$biopsy <-df$biopsy[m]
se$translocation <-df$translocation[m]
se$isotype <-df$isotype[m]
se$site <-df$site[m]
se$source <-df$source[m]

dds <- DESeqDataSet(se, ~1)
ntd <- normTransform(dds)

# Toss "low quality" (as described in original publication) 
dds.discovery <- dds[, dds$cohort == "Discovery"]
ntd.discovery <- normTransform(dds.discovery)

# Write percent variance to tsv
percentVar <- data.frame(attr(pcaData, "percentVar") * 100)
rownames(percentVar)<-paste0('PC',rownames(percentVar))
colnames(percentVar) <- "percent_variance"
write.table(percentVar, "/Users/nate/Projects/EBV_interactome/pca/percent_variance_bl_discovery.tsv", sep="\t")


# # BL Discovery PC3 and PC4 sex
# pcaData <- plotPCA(ntd.discovery, intgroup=c("sex"), returnData=TRUE, pca1=3, pca2=4)
# aa <- ggplot(pcaData) + aes(PC1,PC2,fill=group) + geom_point(size=8,alpha=.7,shape=21)  + theme_bw()
# ggsave('/Users/nate/Projects/EBV_interactome/pca/bl_discovery_sex_pca3_4.svg')
# pcaData <- plotPCA(ntd.discovery, intgroup=c("site"), returnData=TRUE, pca1=3, pca2=4)
# aa <- ggplot(pcaData) + aes(PC1,PC2,fill=group) + geom_point(size=8,alpha=.7,shape=21)  + theme_bw()
# ggsave('/Users/nate/Projects/EBV_interactome/pca/bl_sex_pca3_4.svg')

a <- PCAdf(ntd.discovery, c("ebvstatus"))
pca <- data.frame(a$x)
pca$ebvstatus <- as.numeric(ntd.discovery$ebvstatus) - 1
pca$sex <- as.numeric(ntd.discovery$sex) - 1
pca$genome <- as.numeric(ntd.discovery$genome) - 1
pca$age <- as.numeric(ntd.discovery$age) - 1
pca$biopsy <- as.numeric(ntd.discovery$biopsy) - 1
pca$translocation <- as.numeric(ntd.discovery$translocation) - 1
pca$isotype <- as.numeric(ntd.discovery$isotype) - 1
pca$site <- as.numeric(ntd.discovery$site) - 1
pca$source <- as.numeric(ntd.discovery$source) - 1
pca$variant <- as.numeric(ntd.discovery$variant) - 1

columns <- c("ebvstatus", "sex","age", "genome","biopsy", "translocation", "isotype", "site", "source", "variant" )
rows <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
top_n_PCs <- 10
dd <- data.frame(matrix(nrow=length(rows), ncol=length(columns)))
colnames(dd) <- columns
rownames(dd) <- rows

for (effect in columns) {
  
  li <- c('a','b','c','d','e','f','g','h','i','j') # placeholder list
  ind = 1
  nlevels <- length(unique(pca[[effect]]))
  if (nlevels == 2) {fam='binomial'} else {fam='poisson'}
  print(effect)
  for (pc in c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")){
     forml<-as.formula(paste0(effect, ' ~ ', pc))
     logit <- glm(forml, data=pca, family=fam)
     pred <- predict(logit)
     variant_rsquared <- PseudoR2(logit) # from desctools -> McFadden pseudoR2
     li[ind] <- as.numeric(variant_rsquared)
    #  print(as.numeric(variant_rsquared))
    #  rng = range(pca[[pc]])
    #  span <- rng[2] - rng[1] 
    #  left = rng[2] - 1.5 * span
    #  right = rng[1] + 1.5 * span
    #  xweight <- seq(left, right, .1)
    #  newdata = data.frame('ph' = xweight)
    #  colnames(newdata) <- pc
    #  yweight <- predict(logit, newdata, type="response")
     #pp<-ggplot(pca, aes_string(x=pc, y=effect)) + geom_point(color='black',fill='#ee5622',size=5,alpha=.7,shape=21) + geom_smooth(method="glm", colour='black',method.args=list(family="binomial"),se=FALSE)  + theme_classic() #+  scale_y_continuous(breaks=seq(nlevels) - 1,labels=levels(ntd.discovery[[effect]])) + theme(axis.text.y = element_text(size=16), axis.text.x=element_text(size=16), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) 
     #fh <- paste('/Users/nate/Projects/EBV_interactome/pca/',pc,'_',effect,'_BL.svg', sep='')
     #ggsave(fh)
     ind = ind + 1

  }
  dd[effect] <- li
}

write.table(dd, '/Users/nate/Projects/EBV_interactome/pca/rsquared_pca_logisticreg_BL.tsv', sep='\t')


# Toss formalin fixed samples too
dds.discovery <- dds[, dds$cohort == "Discovery" & dds$biopsy == "FF"]
ntd.discovery <- normTransform(dds.discovery)

# Figure 1A
# BL Discovery PC1 and PC2 EBV status
pcaData <- plotPCA(ntd.discovery, intgroup=c("ebvstatus"), returnData=TRUE, pca1=1, pca2=2) 
percentVar <- round(100 * attr(pcaData, "percentVar"))
aa <- ggplot(pcaData) + aes(PC1,PC2,fill=group) + geom_point(size=8,alpha=.7,shape=21)  + theme_bw() + scale_fill_brewer(palette="Dark2")
ggsave('/Users/nate/Projects/EBV_interactome/pca/bl_noebv_ebvstatus_pca1_2.svg')


# pcaData <- plotPCA(ntd.discovery, intgroup=c("sex"), returnData=TRUE, pca1=3, pca2=4)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# aa <- ggplot(pcaData) + aes(PC1,PC2,fill=group) + geom_point(size=8,alpha=.7,shape=21)  + theme_bw()
# ggsave('/Users/nate/Projects/EBV_interactome/pca/bl_noebv_sex_pca3_4.svg')


a <- plotPCA(ntd.discovery, c("ebvstatus"))
pca <- data.frame(a$x)
pca$ebvstatus <- as.numeric(ntd.discovery$ebvstatus) - 1
pca$sex <- as.numeric(ntd.discovery$sex) - 1
pca$variant <- as.numeric(ntd.discovery$variant) - 1
pca$genome <- as.numeric(ntd.discovery$genome) - 1
pca$age <- as.numeric(ntd.discovery$age) - 1
pca$biopsy <- as.numeric(ntd.discovery$biopsy) - 1
pca$translocation <- as.numeric(ntd.discovery$translocation) - 1
pca$isotype <- as.numeric(ntd.discovery$isotype) - 1
pca$site <- as.numeric(ntd.discovery$site) - 1
pca$source <- as.numeric(ntd.discovery$source) - 1
pca$variant <- as.numeric(ntd.discovery$variant) - 1


dd <- data.frame(matrix(nrow=10,ncol=10))
colnames(dd) <- c("ebvstatus", "sex","age", "genome","biopsy", "translocation", "isotype", "site", "source", "variant" )

for (effect in c("ebvstatus", "sex","age", "genome","biopsy", "translocation", "isotype", "site", "source", "variant" )) {
  
  li <- c('a','b','c','d','e','f','g','h','i','j') # placeholder list
  ind = 1
  nlevels <- length(unique(pca[[effect]]))
  if (nlevels == 2) {fam='binomial'} else {fam='poisson'}
  print(effect)
  for (pc in c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")){
     forml<-as.formula(paste0(effect, ' ~ ', pc))
     logit <- glm(forml, data=pca, family=fam)
     pred <- predict(logit)
     variant_rsquared <- PseudoR2(logit) # from desctools -> McFadden pseudoR2
     li[ind] <- as.numeric(variant_rsquared)
    #  print(as.numeric(variant_rsquared))
    #  rng = range(pca[[pc]])
    #  span <- rng[2] - rng[1] 
    #  left = rng[2] - 1.5 * span
    #  right = rng[1] + 1.5 * span
    #  xweight <- seq(left, right, .1)
    #  newdata = data.frame('ph' = xweight)
    #  colnames(newdata) <- pc
    #  yweight <- predict(logit, newdata, type="response")
    #  pp<-ggplot(pca, aes_string(x=pc, y=effect)) + geom_point(color='black',fill='#ee5622',size=5,alpha=.7,shape=21) + geom_smooth(method="glm", colour='black',method.args=list(family="binomial"),se=FALSE)  + theme_classic() #+  scale_y_continuous(breaks=seq(nlevels) - 1,labels=levels(ntd.discovery[[effect]])) + theme(axis.text.y = element_text(size=16), axis.text.x=element_text(size=16), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) 
    #  fh <- paste('/Users/nate/Projects/EBV_interactome/pca/',pc,'_',effect,'_BL_noebv.svg', sep='')
    #  ggsave(fh)
     ind = ind + 1

  }
  dd[effect] <- li
}
rownames(dd) <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
write.table(dd, '/Users/nate/Projects/EBV_interactome/pca/rsquared_pca_logisticreg_BL_no_FF.tsv', sep='\t')



dds.discovery <- dds[, dds$cohort == "Discovery" & dds$ebvstatus == 'EBV-positive' & dds$biopsy == 'FF']
ntd.discovery <- normTransform(dds.discovery)

a <- plotPCA(ntd.discovery, c("ebvstatus"))
pca <- data.frame(a$x)

pca$ebvstatus <- as.numeric(ntd.discovery$ebvstatus) - 1
pca$sex <- as.numeric(ntd.discovery$sex) - 1
pca$variant <- as.numeric(ntd.discovery$variant) - 1
pca$genome <- as.numeric(ntd.discovery$genome) - 1
pca$age <- as.numeric(ntd.discovery$age) - 1
pca$biopsy <- as.numeric(ntd.discovery$biopsy) - 1
pca$translocation <- as.numeric(ntd.discovery$translocation) - 1
pca$isotype <- as.numeric(ntd.discovery$isotype) - 1
pca$site <- as.numeric(ntd.discovery$site) - 1
pca$source <- as.numeric(ntd.discovery$source) - 1
pca$variant <- as.numeric(ntd.discovery$variant) - 1


dd <- data.frame(matrix(nrow=10,ncol=10))
colnames(dd) <- c("ebvstatus", "sex","age", "genome","biopsy", "translocation", "isotype", "site", "source", "variant" )

for (effect in c("ebvstatus", "sex","age", "genome","biopsy", "translocation", "isotype", "site", "source", "variant" )) {
  
  li <- c('a','b','c','d','e','f','g','h','i','j') # placeholder list
  ind = 1
  nlevels <- length(unique(pca[[effect]]))
  if (nlevels == 2) {fam='binomial'} else {fam='poisson'}
  print(effect)
  for (pc in c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")){
     forml<-as.formula(paste0(effect, ' ~ ', pc))
     logit <- glm(forml, data=pca, family=fam)
     pred <- predict(logit)
     variant_rsquared <- PseudoR2(logit) # from desctools -> McFadden pseudoR2
     li[ind] <- as.numeric(variant_rsquared)
    #  print(as.numeric(variant_rsquared))
    #  rng = range(pca[[pc]])
    #  span <- rng[2] - rng[1] 
    #  left = rng[2] - 1.5 * span
    #  right = rng[1] + 1.5 * span
    #  xweight <- seq(left, right, .1)
    #  newdata = data.frame('ph' = xweight)
    #  colnames(newdata) <- pc
    #  yweight <- predict(logit, newdata, type="response")
    #  pp<-ggplot(pca, aes_string(x=pc, y=effect)) + geom_point(color='black',fill='#ee5622',size=5,alpha=.7,shape=21) + geom_smooth(method="glm", colour='black',method.args=list(family="binomial"),se=FALSE)  + theme_classic() #+  scale_y_continuous(breaks=seq(nlevels) - 1,labels=levels(ntd.discovery[[effect]])) + theme(axis.text.y = element_text(size=16), axis.text.x=element_text(size=16), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) 
    #  fh <- paste('/Users/nate/Projects/EBV_interactome/pca/',pc,'_',effect,'_BL_only_ebvpos.svg', sep='')
    #  ggsave(fh)
     ind = ind + 1

  }
  dd[effect] <- li
}
rownames(dd) <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

write.table(dd, '/Users/nate/Projects/EBV_interactome/pca/bl_only_ebvpos_and_ff.tsv', sep='\t')
