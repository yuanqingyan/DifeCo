#' gbm_dat
#'
#' Genomic alterations data consisting of 23 genes from 232 IDH-WT GBM patients from UThealth cohort and 266 TCGA samples.
#'
#' @usage data(gbm_dat)
#'
#'
#' @keywords datasets
#'
#' @references Yan et al. (2020) JCO Precision Oncology 4, 575-584
#' (\href{https://ascopubs.org/doi/abs/10.1200/PO.19.00385})
#'
#' @examples
#' data(gbm_dat)
#' cohort<-gbm_dat$cohort;unique(cohort)
#' dat<-gbm_dat[,-1]
#'
#'
data(gbm_dat)
