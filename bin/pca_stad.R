library(SummarizedExperiment)
library(readr)
library(DESeq2)
library(ggplot2)
library(vsn)
library(DescTools)



plotPCA = function(object, intgroup="condition", ntop=500, returnData=FALSE,pca1=1,pca2=2)
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


#STOMACH
pheno <- read_tsv('/Users/nate/Projects/EBV_interactome/stad/mmc2.csv')
df <- data.frame(pheno)
df <- data.frame(lapply(df, factor))

counts <-read.csv('/Users/nate/Projects/EBV_interactome/stad_no_ebvgenes.tsv', sep='\t')
counts_matrix = counts[,2:347]
rownames(counts_matrix) <- counts[,1]
counts <- as.matrix(round(counts_matrix,0))
colnames(counts) <- lapply(colnames(counts), function(y) {z=gsub("[.]","-",y);substr(z,1,12)})

se <- SummarizedExperiment(assays=counts)
m <- match(colnames(se), df$TCGA.Participant.Barcode)

se$ebvstatus <- c(rep('EBV negative',320),rep('EBV positive',26))
se$subtype <- df$Molecular_Subtype[m]
se$leukocyte <- df$Leukocyte.Fraction[m]
se$cdkn1a <- df$CDKN2A.methylation[m]
se$gender <-df$Gender[m]
se$race <- df$Race[m]
se$region <- df$Geographic.Region[m]
se$classification <- df$Gastric.histological.classification[m]
se$stage <- df$Pathologic.Stage[m]
se$age <- df$Age.at.initial.pathologic.diagnosis[m]
se$site <- df$Detailed.Anatomic.Site[m]
se$msi <- df$MSI.Status[m]

cspm <- colSums(assay(se))/1e6

dds <- DESeqDataSet(se, ~1)
ntd <- normTransform(dds)

vsd <- vst(dds)
pcaData<-plotPCA(vsd, "ebvstatus",2250,returnData=TRUE, pca1=3, pca2=8) 
aa <- ggplot(pcaData) + aes(PC1,PC2,fill=group) + geom_point(size=8,alpha=.7,shape=21)  + theme_bw() + scale_fill_brewer(palette="Dark2")
ggsave('/Users/nate/Projects/EBV_interactome/pca/ebvstatus_pca3_8_stad.svg')




plotPCA = function(object, intgroup="condition", ntop=500, returnData=FALSE,pca1=1,pca2=2)
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


a <- plotPCA(vsd, c("gender"),2250)
pca <- data.frame(a$x)

pca$ebvstatus <- as.numeric(factor(vsd$ebvstatus)) - 1
logit <- glm(ebvstatus ~ PC3 + PC8, data=pca, family='binomial')
null <- glm(ebvstatus ~ 1, data=pca, family='binomial')
ebv_rsquared <- 1-logLik(logit)/logLik(null)

pca$subtype <- as.numeric(factor(vsd$subtype)) - 1
pca$gender <- as.numeric(factor(vsd$gender)) - 1
pca$site <- as.numeric(factor(vsd$site)) - 1
pca$leukocyte <- as.numeric(factor(vsd$leukocyte)) - 1
pca$stage <- as.numeric(factor(vsd$stage)) - 1
pca$age <- as.numeric(factor(vsd$age))
pca$cdkn1a <- as.numeric(factor(se$cdkn1a))
pca$race <- as.numeric(factor(se$race))  - 1
pca$region <- as.numeric(factor(se$region)) - 1
pca$classification <- as.numeric(factor(se$classification)) - 1
pca$msi <- as.numeric(factor(se$msi)) - 1

dd <- data.frame(matrix(nrow=24,ncol=11))
colnames(dd) <- c("ebvstatus", "gender","age", "stage","classification", "region",  "race", "leukocyte", "site", "subtype", "msi" )

for (effect in colnames(dd)) {
  
  li <- seq(24)#c('a','b','c','d','e','f','g','h','i','j', 'k','l','m','n','o') # placeholder list
  ind = 1
  nlevels <- length(unique(pca[[effect]]))
  if (nlevels == 2) {fam='binomial'} else {fam='poisson'}
  print(effect)
  rown = sprintf("PC%s",seq(24))#("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15")
  for (pc in rown){
     forml<-as.formula(paste0(effect, ' ~ ', pc))
     logit <- glm(forml, data=pca, family=fam)
     pred <- predict(logit)
     variant_rsquared <- PseudoR2(logit) # from desctools -> McFadden pseudoR2
     li[ind] <- as.numeric(variant_rsquared)
     print(as.numeric(variant_rsquared))
     rng = range(pca[[pc]])
     span <- rng[2] - rng[1] 
     left = rng[2] - 1.5 * span
     right = rng[1] + 1.5 * span
     xweight <- seq(left, right, .1)
     newdata = data.frame('ph' = xweight)
     colnames(newdata) <- pc
     yweight <- predict(logit, newdata, type="response")
     pp<-ggplot(pca, aes_string(x=pc, y=effect)) + geom_point(color='black',fill='#ee5622',size=5,alpha=.7,shape=21) + geom_smooth(method="glm", colour='black',method.args=list(family="binomial"),se=FALSE)  + theme_classic() #+  scale_y_continuous(breaks=seq(nlevels) - 1,labels=levels(ntd.discovery[[effect]])) + theme(axis.text.y = element_text(size=16), axis.text.x=element_text(size=16), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) 
     fh <- paste('/Users/nate/Projects/EBV_interactome/pca/',pc,'_',effect,'mssonly_ebvpos_only_STAD.svg', sep='')
     ggsave(fh)
     ind = ind + 1

  }
  dd[effect] <- li
}
rownames(dd) <- rown

write.table(dd, '/Users/nate/Projects/EBV_interactome/pca/stad_pca_r2.tsv', sep='\t')


plotPCA = function(object, intgroup="condition", ntop=500, returnData=FALSE,pca1=1,pca2=2)
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

## Only MSS patients
se <- subset(se, select = colData(se)$msi=="MSS")
dds <- DESeqDataSet(se, ~1)
ntd <- normTransform(dds)
vsd <- vst(dds)
pcaData<-plotPCA(vsd, "ebvstatus",2250,returnData=TRUE, pca1=3, pca2=1) 
aa <- ggplot(pcaData) + aes(PC1,PC2,fill=group) + geom_point(size=8,alpha=.7,shape=21)  + theme_bw() + scale_fill_brewer(palette="Dark2")
ggsave('/Users/nate/Projects/EBV_interactome/pca/ebvstatus_pca3_1_stad.svg')


plotPCA = function(object, intgroup="condition", ntop=500, returnData=FALSE,pca1=1,pca2=2)
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


a <- plotPCA(vsd, c("gender"),2250)
pca <- data.frame(a$x)

pca$ebvstatus <- as.numeric(factor(vsd$ebvstatus)) - 1
logit <- glm(ebvstatus ~ PC3 + PC8, data=pca, family='binomial')
null <- glm(ebvstatus ~ 1, data=pca, family='binomial')
ebv_rsquared <- 1-logLik(logit)/logLik(null)

pca$subtype <- as.numeric(factor(vsd$subtype)) - 1
pca$gender <- as.numeric(factor(vsd$gender)) - 1
pca$site <- as.numeric(factor(vsd$site)) - 1
pca$leukocyte <- as.numeric(factor(vsd$leukocyte)) - 1
pca$stage <- as.numeric(factor(vsd$stage)) - 1
pca$age <- as.numeric(factor(vsd$age))
pca$cdkn1a <- as.numeric(factor(se$cdkn1a))
pca$race <- as.numeric(factor(se$race))  - 1
pca$region <- as.numeric(factor(se$region)) - 1
pca$classification <- as.numeric(factor(se$classification)) - 1
pca$msi <- as.numeric(factor(se$msi)) - 1

dd <- data.frame(matrix(nrow=25,ncol=12))
colnames(dd) <- c("ebvstatus", "gender","age", "stage","classification", "region",  "cdkn1a","race", "leukocyte", "site", "subtype", "msi" )

for (effect in colnames(dd)) {
  
  li <- seq(25)#c('a','b','c','d','e','f','g','h','i','j', 'k','l','m','n','o') # placeholder list
  ind = 1
  nlevels <- length(unique(pca[[effect]]))
  if (nlevels == 2) {fam='binomial'} else {fam='poisson'}
  print(effect)
  rown = sprintf("PC%s",seq(25))#("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15")
  for (pc in rown){
     forml<-as.formula(paste0(effect, ' ~ ', pc))
     logit <- glm(forml, data=pca, family=fam)
     pred <- predict(logit)
     variant_rsquared <- PseudoR2(logit) # from desctools -> McFadden pseudoR2
     li[ind] <- as.numeric(variant_rsquared)
     print(as.numeric(variant_rsquared))
     rng = range(pca[[pc]])
     span <- rng[2] - rng[1] 
     left = rng[2] - 1.5 * span
     right = rng[1] + 1.5 * span
     xweight <- seq(left, right, .1)
     newdata = data.frame('ph' = xweight)
     colnames(newdata) <- pc
     yweight <- predict(logit, newdata, type="response")
     pp<-ggplot(pca, aes_string(x=pc, y=effect)) + geom_point(color='black',fill='#ee5622',size=5,alpha=.7,shape=21) + geom_smooth(method="glm", colour='black',method.args=list(family="binomial"),se=FALSE)  + theme_classic() #+  scale_y_continuous(breaks=seq(nlevels) - 1,labels=levels(ntd.discovery[[effect]])) + theme(axis.text.y = element_text(size=16), axis.text.x=element_text(size=16), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) 
     fh <- paste('/Users/nate/Projects/EBV_interactome/pca/',pc,'_',effect,'mssonly_STAD.svg', sep='')
     ggsave(fh)
     ind = ind + 1

  }
  dd[effect] <- li
}
rownames(dd) <- rown

write.table(dd, '/Users/nate/Projects/EBV_interactome/pca/stad_mss_only.tsv', sep='\t')




plotPCA = function(object, intgroup="condition", ntop=500, returnData=FALSE,pca1=1,pca2=2)
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

## Only MSS patients
se <- subset(se, select = c(colData(se)$msi=="MSS" & colData(se)$ebvstatus == "EBV Positive"))
dds <- DESeqDataSet(se, ~1)
ntd <- normTransform(dds)
vsd <- vst(dds)
pcaData<-plotPCA(vsd, "ebvstatus",2250,returnData=TRUE, pca1=3, pca2=1) 
aa <- ggplot(pcaData) + aes(PC1,PC2,fill=group) + geom_point(size=8,alpha=.7,shape=21)  + theme_bw() + scale_fill_brewer(palette="Dark2")
ggsave('/Users/nate/Projects/EBV_interactome/pca/ebvstatus_pca3_1_stad.svg')


plotPCA = function(object, intgroup="condition", ntop=500, returnData=FALSE,pca1=1,pca2=2)
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


a <- plotPCA(vsd, c("gender"),2250)
pca <- data.frame(a$x)

pca$ebvstatus <- as.numeric(factor(vsd$ebvstatus)) - 1
logit <- glm(ebvstatus ~ PC3 + PC8, data=pca, family='binomial')
null <- glm(ebvstatus ~ 1, data=pca, family='binomial')
ebv_rsquared <- 1-logLik(logit)/logLik(null)

pca$subtype <- as.numeric(factor(vsd$subtype)) - 1
pca$gender <- as.numeric(factor(vsd$gender)) - 1
pca$site <- as.numeric(factor(vsd$site)) - 1
pca$leukocyte <- as.numeric(factor(vsd$leukocyte)) - 1
pca$stage <- as.numeric(factor(vsd$stage)) - 1
pca$age <- as.numeric(factor(vsd$age))
pca$cdkn1a <- as.numeric(factor(se$cdkn1a))
pca$race <- as.numeric(factor(se$race))  - 1
pca$region <- as.numeric(factor(se$region)) - 1
pca$classification <- as.numeric(factor(se$classification)) - 1
pca$msi <- as.numeric(factor(se$msi)) - 1

dd <- data.frame(matrix(nrow=25,ncol=12))
colnames(dd) <- c("ebvstatus", "gender","age", "stage","classification", "region",  "cdkn1a","race", "leukocyte", "site", "subtype", "msi" )

for (effect in colnames(dd)) {
  
  li <- seq(25)#c('a','b','c','d','e','f','g','h','i','j', 'k','l','m','n','o') # placeholder list
  ind = 1
  nlevels <- length(unique(pca[[effect]]))
  if (nlevels == 2) {fam='binomial'} else {fam='poisson'}
  print(effect)
  rown = sprintf("PC%s",seq(25))#("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15")
  for (pc in rown){
     forml<-as.formula(paste0(effect, ' ~ ', pc))
     logit <- glm(forml, data=pca, family=fam)
     pred <- predict(logit)
     variant_rsquared <- PseudoR2(logit) # from desctools -> McFadden pseudoR2
     li[ind] <- as.numeric(variant_rsquared)
     print(as.numeric(variant_rsquared))
     rng = range(pca[[pc]])
     span <- rng[2] - rng[1] 
     left = rng[2] - 1.5 * span
     right = rng[1] + 1.5 * span
     xweight <- seq(left, right, .1)
     newdata = data.frame('ph' = xweight)
     colnames(newdata) <- pc
     yweight <- predict(logit, newdata, type="response")
     pp<-ggplot(pca, aes_string(x=pc, y=effect)) + geom_point(color='black',fill='#ee5622',size=5,alpha=.7,shape=21) + geom_smooth(method="glm", colour='black',method.args=list(family="binomial"),se=FALSE)  + theme_classic() #+  scale_y_continuous(breaks=seq(nlevels) - 1,labels=levels(ntd.discovery[[effect]])) + theme(axis.text.y = element_text(size=16), axis.text.x=element_text(size=16), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) 
     fh <- paste('/Users/nate/Projects/EBV_interactome/pca/',pc,'_',effect,'mssonly_ebvposonly_STAD.svg', sep='')
     ggsave(fh)
     ind = ind + 1

  }
  dd[effect] <- li
}
rownames(dd) <- rown

write.table(dd, '/Users/nate/Projects/EBV_interactome/pca/stad_mss_only_ebvpos_only.tsv', sep='\t')

