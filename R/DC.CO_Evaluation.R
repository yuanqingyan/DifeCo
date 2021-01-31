#' @title DC.CO_Evaluation
#'
#' @description Evaluate the cooccurrence/mutual exclusion or differential cooccurrence/mutual exclusion of the gene mutation data
#'
#' @param input_data Dataframe or matrix of genomic alteration data (for example: gene mutation data). Should be a binary data with 0 as the wild type while 1 for altered/mutation.
#' @param mode Which mode to be used. Three modes are included: "DC"(evaluate the significance of differential cooccurrence between two different groups),"Separate"(evaluate the pairwise cooccurrence in each group separately),"Single"(evalute the pairwise cooccurrence for the entire dataset).
#' @param group Binary data.  Applys to "Separate" or "DC" mode.
#' @param which_group_to_be_one One category of the group parameter. Applys to "Separate" or "DC" mode. For "DC" mode, this category is modeled.
#' @param adjust.method Method for p value adjustment. Should be one of the "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none".
#' @param FDRCutoff FDR cutoff
#'
#'
#'
#' @return A list with odds ratio, raw p value and fasle discovery rate
#' @export
#' @examples
#' data(gbm_dat)
#' cohort<-gbm_dat$cohort;unique(cohort)
#' dat<-gbm_dat[,-1]
#'
#' #Evaluate pairwise cooccurrence of genomic alteration of the entire gbm_dat data.
#' Result_Single<-DC.CO_Evaluation(input_data=dat,
#'                                 mode="Single",
#'                                 adjust.method="BH",
#'                                 FDRCutoff=0.05)
#'
#' #Evaluate pairwise cooccurrence for the two cohorts separately.
#' Result_Separate<-DC.CO_Evaluation(input_data=dat,
#'                                   group=cohort,
#'                                   which_group_to_be_one="UT",
#'                                   mode="Separate",
#'                                   adjust.method="BH",
#'                                   FDRCutoff=0.05)
#'
#' #Evaluate the significance of differential cooccurrence between patients with PTEN=0 and PTEN=1
#' UT_dat<-gbm_dat[gbm_dat$cohort=="UT",-1]
#' PTEN_grp<-UT_dat$PTEN;table(PTEN_grp)
#' dat2<-subset(UT_dat,select=-PTEN)
#' Result_DC<-DC.CO_Evaluation(input_data=dat2,
#'                            group=PTEN_grp,
#'                            which_group_to_be_one=1,
#'                            mode="DC",
#'                            adjust.method="BH",
#'                            FDRCutoff=0.1)



DC.CO_Evaluation<-function(input_data=input_data,
                           group=NULL,
                           mode="DC",
                           which_group_to_be_one=NULL,
                           adjust.method="BH",
                           FDRCutoff=0.05){
  if(!(is.data.frame(input_data) | is.matrix(input_data))){
    stop("Input data should be in a data frame or matrix!")
  }
  if(!(mode %in% c("DC","Separate","Single"))){
    stop("mode should be one of 'DC','Separate' or 'Single'!")
  }
  if(mode %in% c("DC","Separate")){
    if(!(length(unique(group))==2)){
      stop("Only two groups are allowed!")
    }
    if(!(which_group_to_be_one %in% unique(group))){
      stop("which_group_to_be_one is not in group!")
    }
  }

  if(!adjust.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")){
    stop("adjust.method should be one of ['holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY','fdr', 'none']!")
  }

  if(!is.numeric(FDRCutoff)){
    stop("FDRCutoff should be numeric data and ranges from 0~1!")
  }
  input_data<-input_data[,!apply(input_data,2,function(x) all(is.na(x)))]
  if(mode=="DC"){
    #### calculate inter and main P value OR value
    P_OR.list<-calculate_sig_interaction(mat=input_data,
                                         group=group,
                                         which_group_to_be_one=which_group_to_be_one,
                                         adjust.method=adjust.method)
    ##stpe 2## p value adjustment
    list.P.adj<-adjust_inter_P(inputList=P_OR.list,
                               FDRCutoff=FDRCutoff,
                               adjust.method=adjust.method)
  }
  if(mode=="Separate"){
    #### calculate single group: P value OR value
    P_OR.list.sep.group<-make_P_OR_SeqGroup(mat=input_data,
                                            group=group,
                                            which_group_to_be_one=which_group_to_be_one,
                                            adjust.method=adjust.method)
    ##stpe 2## p value adjustment
    list.P.adj<-adjust_SepGroup_P(inputList=P_OR.list.sep.group,
                                  FDRCutoff=FDRCutoff,
                                  adjust.method=adjust.method)
  }
  if(mode=="Single"){
    #### calculate single group: P value OR value
    P_OR.list.single<-make_P_OR_Single(mat=input_data,adjust.method=adjust.method)
    ##stpe 2## p value adjustment
    list.P.adj<-adjust_Single_P(inputList=P_OR.list.single,
                                FDRCutoff=FDRCutoff,
                                adjust.method=adjust.method)
    group=which_group_to_be_one="Single"
  }
  list.P.adj$group=group
  list.P.adj$which_group_to_be_one=which_group_to_be_one
  list.P.adj$mode=mode
  return(list.P.adj)
}
