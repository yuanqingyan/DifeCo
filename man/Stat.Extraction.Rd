% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stat.Extraction.R
\name{Stat.Extraction}
\alias{Stat.Extraction}
\title{Stat.Extraction}
\usage{
Stat.Extraction(obj = obj_input, WhichGrp_adjP.Separate = NULL)
}
\arguments{
\item{obj}{Object generated by DC.CO_Evaluation function. Required}

\item{WhichGrp_adjP.Separate}{Which group to adjust p value in separate mode. Restricts to separate mode only. Option includes: NULL(p value adjusted by the combined dataset from two groups), "All"(p value adjusted by the combined dataset from two groups), "None"(no p value adjustment), "Each"(p value adjusted within each group). The value can also be one of the category of the group parameter. In this case, p value within the specified group is adjusted while the other group shows the raw p value. Optional. Default: NULL}
}
\value{
A list with multiple objects. The statistics of the result can be obtained from the element with the name of "Stat". The element of "Stat" is a data frame with log odds ratio(logOR), raw p value(RawP), fasle discovery rate(FDR) and significance. LogOR>0 implies the genes are in an cooccurrence pattern, while <0 for mutual exclusivity.
}
\description{
Extract the statistics from the object generated by DC.CO_Evaluation
}
\examples{
data(gbm_dat)
dat<-gbm_dat[,-1]
Result_Single<-DC.CO_Evaluation(input_data=dat,
                                mode="Single",
                                adjust.method="BH",
                                FDRCutoff=0.05)
sta_Sin<-Stat.Extraction(obj=Result_Single)
sig_Sin<-sta_Sin$Stat;head(sig_Sin)
sig_Sin[sig_Sin$Sig=="Yes",]


cohort<-gbm_dat$cohort;unique(cohort)
dat<-gbm_dat[,-1]

Result_Separate<-DC.CO_Evaluation(input_data=dat,
                                  group=cohort,
                                  which_group_to_be_one="UT",
                                  mode="Separate",
                                  adjust.method="BH",
                                  FDRCutoff=0.05)
sta_Sep<-Stat.Extraction(obj=Result_Separate, WhichGrp_adjP.Separate="Each")
sig_Sep<-sta_Sep$Stat;head(sig_Sep)
sig_Sep[sig_Sep$Sig.In.UT=="Yes",]

}
