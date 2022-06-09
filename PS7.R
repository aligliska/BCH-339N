## BCH 339N
## RNA-Seq Analysis - 21pt

## Alice Gee, ag67642

## We learned about NGS and RNA-Seq in class.
## In this assignment, we will analyze RNA-Seq data using DESeq2. 
## If you want to learn more about this package, check out the vignette
## https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
## We are going to carry out exploratory and differential expression analyses. 

library('DESeq2')

## In this homework, we will use the fread function from the data.table package. 
## fread allows us to extract .csv.gz files directly from the web. 
## If you are interested in learning more, data.table is a very popular R package.
library('R.utils')
library('data.table')

## We are going to analyze the following dataset
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65912
## Note that there is a link to .csv.gz file with RNA Expression on this page.
## https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65912/suppl/GSE65912_RNASeq.csv.gz

rnaseq_counts = fread ( 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65912/suppl/GSE65912_RNASeq.csv.gz', 
                        data.table = F)

## We will assign gene names to the row.names attribute below. 
## The column headers correspond to samples and look like: "GM12878_Rep1_Counts"
## GM12878 corresponds to the cell line identifier. 
## Rep1 corresponds to the replicate number
row.names(rnaseq_counts) = rnaseq_counts$V1
rnaseq_counts = as.matrix (rnaseq_counts[,-1] ) 

## Q1 - 6pt 
## Data Exploration: 

## A typical exploratory analysis is to inspect a scatterplot of read counts
## Plot two scatterplots using log10 of the gene read counts. 
## For the first plot compare GM12878_Rep1_Counts to GM12878_Rep2_Counts
## For the second plot compare GM12878_Rep1_Counts to GM12891_Rep1_Counts
## Describe your observations regarding these two plots in a few sentences. 
log_rnaseq_counts = log10(rnaseq_counts + 1)
plot(log_rnaseq_counts[, "GM12878_Rep1_Counts"], log_rnaseq_counts[, "GM12878_Rep2_Counts"])
abline(lm(log_rnaseq_counts[, "GM12878_Rep1_Counts"] ~ log_rnaseq_counts[, "GM12878_Rep2_Counts"]), col = "red")
title("GM12878_Rep1_Counts to GM12878_Rep2_Counts")

plot(log_rnaseq_counts[, "GM12878_Rep1_Counts"], log_rnaseq_counts[, "GM12891_Rep1_Counts"])
abline(lm(log_rnaseq_counts[, "GM12878_Rep1_Counts"] ~ log_rnaseq_counts[, "GM12891_Rep1_Counts"]), col = "blue")
title("GM12878_Rep1_Counts to GM12891_Rep1_Counts")

# From the scatter plot, we can see that there is a general positive trend across GM12878_Rep1_Counts and GM12878_Rep2_Counts, 
# and across GM12878_Rep1_Counts and GM12891_Rep1_Counts. That said, there is less variance in the first scatter plot
# compared to the second scatter plot. For the most part, as the counts increased for GM12878_Rep1_Counts, there was also 
# an increase in counts across GM12878_Rep2_Counts. Similar in the second plot, as the counts for GM12878_Rep1_Counts
# increased, the counts for GM12891_Rep1_Counts also increased, but there are more points that deviate from the best
# fit line (i.e. more variance). 

## Q2 - 3pt
## Calculate all pairwise Spearman correlations between all the samples
## Draw a heatmap of the calculated pairwise Spearman correlations
## Write a few lines describing the pattern(s) you observe in the heatmap
## The key functions that you will need are provided below.
library('pheatmap')

cormat <- rnaseq_counts %>% cor(use = "pair", method = "spearman")
pheatmap(cormat)
# In the heatmap, we can see that most of the Spearman correlations are relatively high, with the lowest correlation
# values being around 0.95. In the center, we see a diagonal line indicating the highest correlation of 1.0, which 
# aligns with the correlation with itself (e.g. a correlation between GM12878_Rep1_Counts and GM12878_Rep1_Counts). 
# As we move away from the diagonal line of highest correlation, the correlation values decrease. Around the diagonal 
# line of correlation = 1, there are yellow boxes, which form a 3x3 box, that represent the correlations of the replicates
# for each cell line sample. Essentially across the heatmap, we can see the gradiation change by 3x3 sections, aligning
# with each cluster/group of cell line samples and the respective replicates. The lowest correlations are situated
# in the bottom left and upper right corners of the heatmap, which are essentially mirrored correlation comparisons
# of each other (e.g. GM19238_Rep3_Counts and GM12892_Rep2_Counts vs. GM12892_Rep2_Counts and GM19238_Rep3_Counts). 
# Overall, we can see that the concentration of highest correlated values are near the diagonal line from the top left
# to the bottom right corner. A gradiation from highest to lower correlation disperses from the diagonal line to the 
# bottom left and top right corners. 


## We will  extract the cell line names from the column names. 
## We will then define a factor with two levels of length equal to the number of columns of rnaseq_counts
##  the reference level of this factor will match the cell lines:
## c("GM19238", "GM19240", "GM12892", "GM12878" )
## The alternative level of the factor will match: c("GM19239", "GM12891" ) 

cell_line_ids = sapply (strsplit(colnames(rnaseq_counts), "_"), "[[" , 1 ) 
Factor_level_1 = c("GM19238", "GM19240", "GM12892", "GM12878" )
Factor_level_2 =  c("GM19239", "GM12891" )
factor_of_interest = cell_line_ids %in% Factor_level_1
factor_of_interest = as.factor(factor_of_interest)
factor_of_interest = relevel(factor_of_interest, ref = "TRUE")


## Use DESeqDataSetFromMatrix to prepare the data for differential expression analysis.
dds <- DESeqDataSetFromMatrix (countData = rnaseq_counts,
                               colData = DataFrame(Variable = factor_of_interest ),
                               design = ~ Variable)

counts_per_million <- function (count_matrix) { 
  
  return (t(apply(count_matrix, 1, function(x){1000000*x/colSums(rnaseq_counts)})))
  
}

## Q3 - 3pt
## Next, we will filter out low expressed genes from our DESeq object
## Specifically, we will only keep the genes that have a counts_per_million greater than 1 in more than 3 samples. 
count_matrix <- counts(dds)
output <- counts_per_million(count_matrix)
keep <- c()
for (i in 1:nrow(output)){
  # sum across the rows for all that are greater than 1
  temp <- sum(output[i, ] > 1)
  # if there are more than 3 samples that are greater than 1
  if (temp > 3){
    keep = c(keep, i)
  } 
}
dds <- dds[keep,]

## Q4 - 3pt
## We will  run DESeq to identify differentially expressed genes
dds <- DESeq(dds)
resultsNames(dds)
results(dds)
resLFC <- lfcShrink(dds, coef= "Variable_FALSE_vs_TRUE")
resLFC <- resLFC[order(resLFC$pvalue),]
summary(resLFC, alpha = 0.05)

## In the above code, what does lfcShrink achieve? Explain in your own words.
# lfcShrink is a function that changes the log2FoldChange and lfcSE by shrinking large fold change values. 
# The function does this by looking at large log2FoldChange values that are caused by low count values and 
# shrinks them such that log fold changes can be compared across experiments and it is easier to make 
# visualizations. Essentially, lfcShrink shrinks large log fold change values by reassessing the overall 
# distribution of the log fold changes, which in terms changes the lfcSE. In this case in particular, 
# we are shrinking over the coefficient "Variable_FALSE_vs_TRUE" to make estimates for the model matrix
# following the distribution of the maximum likelihood estimation method (i.e. more of a functional description). 

## Q5 - 3pt

## Generate an MA-plot using resLFC

# different y-axis scales so you can see the plot from the zoomed in version and the more full distribution.   
plotMA(resLFC)
plotMA(resLFC, ylim=c(-5,5))

## Q6 -3pt
## Next we will create an interactive html to explore our results
library("ReportingTools")
des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
                         title = 'RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./report")
## This might take a while. 
publish(dds,des2Report, pvalueCutoff=0.05,
        annotation.db="org.Hs.eg", factor = factor_of_interest,
        reportDir="./reports")
finish(des2Report)

## If everything went well, you should be able to inspect an html file containing your results.
## If successful, include one boxplot generated from this report. 
## Even if the above code fails, you can see the most differentially expressed genes with the following code
resLFC


## The top five differentially expressed genes should be: 
## KDM5C, PSMA8, KDM6A, ZFX, SMC1A
## Search what each of these genes are. 
## Do you notice any shared feature(s)? If so, can you speculate what our factor_of_interest was? 

# As expected, the top five differentially expressed genes are KDM5C, PSMA8, KDM6A, ZFX, SMC1A. 
# Of these genes, 4/5 (i.e. KDM5C, KDM6A, ZFX, and SMC1A) of these genes have a decrease in expression levels, 
# as indicated by the negative log2FoldChange. These same 4/5 genes had similar log fold change standard error 
# between 0.05-0.068. The PSMA8 gene was different from the other 4 genes in that it is significantly expressed
# more compared to normal levels (i.e. positive log2 fold change) and had a log fold change standard error of 0.389.
# Overall, we should note that the data used to run differential analysis come from lymphoblastoid cells, which are 
# derived from lymphocytes (which are a type of white blood cell). When looking up the role of each of the five
# differentially expressed genes, we see that KDM5C and KDM6A both code for Lysine-specific demethylases, which 
# play pivotal roles in coding histone demethylase enzymes that help regulate gene transcription. PSMA8 codes for 
# proteasomes that promote acetylation-dependent degradation of histones, which in the end also plays a role in 
# regulating gene transcription. ZFX is a protein-coding gene that plays an important role in cell growth, and 
# is also associated with being a form of transcriptional activator. The SMC1A gene regulates the structure and 
# organization of chromosome 1. As a generalization, all of these genes play a significant role in coding regulatory 
# genes or of similar importance. Mutations at these genes would likely cause detrimental effects on a person. 
# Looking at these genes overall, it appears that we are interested in looking at the variation of RNA levels
# across different genotypes, expression of transcription-related genes, and how it in turn affects the overall
# translation of genes (i.e. leading to gene expression). By focusing on lymphoblastoid cells, it appears that 
# we are interested in looking at the variation of gene expression across different genotypes/mutants to 
# better elucidate the importance of personalized medicine. 
