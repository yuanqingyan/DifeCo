#' @title gbm_dat
#' @name gbm_dat
#' @description Genomic alterations data from IDH-WT GBM patients. The dataset consists of 23 genes. A total of 498 samples from 2 different cohort are involved (232 from UThealth at houston and 266 from TCGA cohort).Please refer to the publication for more details.
#' @docType data
#' @usage data(gbm_dat)
#' @keywords datasets
#' @references Yan et al. (2020) JCO Precision Oncology 4, 575-584
#' (\href{JCOPO}{https://ascopubs.org/doi/abs/10.1200/PO.19.00385})
#' @examples
#' data(gbm_dat)
#' cohort<-gbm_dat$cohort
#' dat<-gbm_dat[,-1]

data(gbm_dat)
