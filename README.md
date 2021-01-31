# DifeCo: Differential Cooccurrence or Mutual Exclusion of Binary Genomic Alteration Data

## Description
Gene mutation resulting in functional dysregulation is the direct cause of most genetic diseases. In many diseases, some gene mutations are tend to occur together and compensate the biological functions with each other. While for some other mutations, their functions are redundant and tend to mutually exclude with each other. This phenomenon is commonly seen in cancer biology. For example, in IDH-WT GBM, TP53 and RB1 mutations often cooccurring, while CDKN2A/B loss is mutually exclusive with TP53 mutation. Due to the disease heterogeneousness, the pattern of gene mutation cooccurrence/mutual exclusion could vary, such as patients with long vs short survival, or patients between different subtypes. The differential cooccurrence/mutual exclusion of gene mutations could be critical for disease treatment. DifeCo is an R package to evaluate the differential occurrence/mutual exclusion of gene mutation. It fits a Firth's bias-reduced logistic regression model between pairwise genes plus the additional group variable. An interaction term of independent predictors is introduced and its significance is evaluated. After the multiplicity adjustment, the pairs of gene are regarded to be statistically significant if the adjusted p value of interaction term is less than the designed cutoff. For the model with interaction term failing to reach significance, the additive model without interaction term is fit to evaluate the cooccurrence/mutual exclusion in the entire dataset. In addition to test the differential cooccurrence/mutual exclusion, DifeCo package can also be used to evaluate and visualize pairwise gene cooccurrence/mutual exclusion in two datasets (Separate mode) or single one dataset (Single mode). Which model to be used purely depends on the hypothesis as well as the nature of the data.
 

## Installation
Installing DifeCo from GitHub
```{r, eval=FALSE}
library(devtools)
install_github("yuanqingyan/DifeCo")
```
## Load data
Load IDH-WT gbm data. Here we focus on UT cohort and assume the patients can be splitted into two groups based on the status of PTEN. Here we have PTEN-WT and PTEN-Alt. This example is to illustrate how to use package with DC mode. 
```{r,eval=TRUE}
library(DifeCo)
data(gbm_dat)
head(gbm_dat[,1:5])
MutDat_UT<-gbm_dat[gbm_dat$cohort=="UT",]
#Remove cohort clumn
MutDat_UT<-MutDat_UT[,-1]
PTEN<-MutDat_UT$PTEN;table(PTEN)
MutDat_WoPTEN<-subset(MutDat_UT,select=-PTEN);head(MutDat_WoPTEN[,1:5])

```
## DC mode
First step is to evaluate the differential cooccurrence/mutual exclusion of genomic alterations between the patients with wild type PTEN (PTEN-WT) and PTEN function altered (PTEN-Alt). We set up the FDR cutoff equal to 0.1. 
```{r, eval=TRUE}
  Result_DC<-DC.CO_Evaluation(input_data=MutDat_WoPTEN,
                              group=PTEN,
                              which_group_to_be_one=1,
                              mode="DC",
                              adjust.method="BH",
                              FDRCutoff=0.1)
```
The gene pairs with significantly differential cooccurrence/mutual exclusion can be extracted by Stat.Extraction function. In this study, the pair of PIK3CA and PIK3R1 shows significantly differential cooccurrence/mutual exclusion. In other word, the pattern of cooccurrence/mutual exclusion of PIK3CA and PIK3R1 depends on the group. In details, alterations in PIK3CA and PIK3RI are mutually excluded in PTEN-WT group, while they cooccur in PTEN-Alt group. For the gene pairs without significant differential cooccurrence/mutual exclusion, the pattern is evaluated in the entire dataset and 13 pairs show significant.
```{r, eval=TRUE}
  sta_DC<-Stat.Extraction(obj=Result_DC)
  sig_DC<-sta_DC$Stat
#Extract the gene pairs with significant differential cooccurrence/mutual exclusion
  sig_DC[sig_DC$Sig.In.DC=="Yes",]
#Extract the gene pairs with significant cooccurrence/mutual exclusion in UT chort dataset
  nrow(sig_DC[sig_DC$Sig.In.CO=="Yes",])
  head(sig_DC[sig_DC$Sig.In.CO=="Yes",])
```
The following code is to visualize the result. The grid with green color shows the gene pair with significant differential cooccurrence/mutual exclusion. The cells with * are the ones with significant coocurrence/mutual exclusion for entire UT dataset.
```{r, eval=TRUE}
  DC.CO_plot(obj=sta_DC,
             label.gene.cex=0.8)
```

## Separate mode
Significance of cooccurrence/mutual exclusion is evaluated in each sub dataset in "Separate" mode. This is particularly helpful to investigate the gene pair patterns in two different datasets (One dataset for discovery purpose, while the other for validation). The example below evaluates cooccurrence/mutual exclusion in UT cohort and validate the result in TCGA cohort. 
```{r, eval=TRUE}
  Cohort=gbm_dat$cohort
  MutDat<-gbm_dat[,-1]
  Result_Sep<-DC.CO_Evaluation(input_data=MutDat,
                               group=Cohort,
                               which_group_to_be_one='UT',
                               mode="Separate",
                               adjust.method="BH",
                               FDRCutoff=0.05)
```
Extract the significant gene pairs. Here, we use WhichGrp_adjP.Separate="Each" to do the multiplicity adjustment in UT and TCGA dataset separately.
```{r, eval=TRUE}
  sta_Sep<-Stat.Extraction(obj=Result_Sep,
                           WhichGrp_adjP.Separate="Each")
  sig_Sep<-sta_Sep$Stat
#Extract the significant gene pairs in UT cohort
  sig_Sep[sig_Sep$Sig.In.UT=="Yes",]
```
The result can be plotted as below.
```{r, eval=TRUE}
  DC.CO_plot(obj=sta_Sep,
             label.gene.cex=0.8,
             sig_CO.cex = 1.5)
```


## Single mode
Assuming data from UT and TCGA cohort can be simply combined into one dataset, we analyze the cooccurrence and mutual exclusion of the combined dataset in "Single" mode.
```{r, eval=TRUE}
 #Remove "cohort" column
  CombineDat<-gbm_dat[,-1]
  Result_Sin<-DC.CO_Evaluation(input_data=CombineDat,
                               mode="Single",
                               adjust.method="BH",
                               FDRCutoff=0.05)
```
Extract the significant gene pairs.
```{r, eval=TRUE}
  sta_Sin<-Stat.Extraction(obj=Result_Sin)
  sig_Sin<-sta_Sin$Stat
#Extract the significant gene pairs
  sig_Sin[sig_Sin$Sig=="Yes",]
```
The result is plotted.
```{r, eval=TRUE}
  DC.CO_plot(obj=sta_Sin,
             label.gene.cex=0.8)
```
