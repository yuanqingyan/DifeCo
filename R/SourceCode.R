
extract_stat_single<-function(obj=obj){
  temp<-obj$CoordList$Single
  temp_ElogOR<-obj$ElogOR
  colN_temp<-colnames(temp_ElogOR)
  rowN_temp<-rownames(temp_ElogOR)

  Col_index<-temp$X
  Row_index<-nrow(temp_ElogOR)+1-temp$Y

  out<-data.frame(Gene1=rowN_temp[Row_index],
                  Gene2=colN_temp[Col_index],
                  temp[,3:5])
  colnames(out)[3:5]<-c("LogOR","Raw.P","FDR")
  row.names(out)<-1:nrow(out)
  out$Sig<-ifelse(out$FDR<obj$FDRCutoff,"Yes","No")
  return(out)
}


extract_stat_DC<-function(obj=Result_DC){
  uniGroup<-unique(obj$group)
  which_group_to_be_one<-obj$which_group_to_be_one
  name_grp0<-uniGroup[!uniGroup==which_group_to_be_one]
  name_grp1<-which_group_to_be_one

  FDR_cutoff<-obj$FDRCutoff
  temp<-get_coordList(obj=obj,n=7)
  pval_DE<-temp$pvalue_DE
  FDR_DE<-temp$FDR_DE
  OR_0<-temp$OR_grp0
  OR_1<-temp$OR_grp1
  OR_CO<-temp$OR_CO
  pval_CO<-temp$pvalue_CO
  FDR_CO<-temp$FDR_CO

  temp_OR<-obj$OR_CO
  colN_temp<-colnames(temp_OR)
  rowN_temp<-rownames(temp_OR)

  Col_index<-OR_CO$X
  Row_index<-nrow(temp_OR)+1-OR_CO$Y

  out<-data.frame(Gene1=rowN_temp[Row_index],
                  Gene2=colN_temp[Col_index],
                  LogOR_Group0=OR_0[,3],
                  LogOR_Group1=OR_1[,3],
                  RawP_Test.DC=pval_DE[,3],
                  FDR_Test.DC=FDR_DE[,3],
                  LogOR_AdjustGroup=OR_CO[,3],
                  RawP_AdjustGroup=OR_CO[,3],
                  FDR_AdjustGroup=FDR_CO[,3])
  row.names(out)<-1:nrow(out)

  for(i in 1:nrow(out)){
    #############
    if(is.na(out$FDR_Test.DC[i])){
      out$LogOR_Group0[i]<-out$LogOR_Group1[i]<-out$RawP_Test.DC[i]<-out$FDR_Test.DC[i]<-NA
    }else{
      if(out$FDR_Test.DC[i]>=FDR_cutoff){out$LogOR_Group0[i]<-out$LogOR_Group1[i]<-out$RawP_Test.DC[i]<-out$FDR_Test.DC[i]<-NA}
    }
    ###############
    if(!is.na(out$FDR_Test.DC[i])){
      out$LogOR_AdjustGroup[i]<-out$RawP_AdjustGroup[i]<-out$FDR_AdjustGroup[i]<-NA
    }
  }
  colnames(out)[3:4]<-sprintf("LogOR_WithinGroup:%s",c(name_grp0,name_grp1))
  out$Sig.In.DC<-ifelse(!is.na(out$FDR_Test.DC),"Yes","No")
  out$Sig.In.CO<-ifelse(!is.na(out$FDR_AdjustGroup) & out$FDR_AdjustGroup<FDR_cutoff,"Yes","No")

  return(out)
}




extract_stat_Sep<-function(obj=Result_Separate,WhichGrp_adjP.Separate=NULL){
  obj_coord<-make_Coord_df(obj=obj,
                           DE_mode=obj$mode,
                           FDRCutoff=obj$FDRCutoff,
                           col_negOR="blue",
                           col_zero="white",
                           col_posOR="red")

  uniGroup<-unique(obj$group)
  which_group_to_be_one<-obj$which_group_to_be_one
  name_grp0<-uniGroup[!uniGroup==which_group_to_be_one]
  name_grp1<-which_group_to_be_one

  sep_checkSig<-check_Sig_SepGrp(dataIn=obj_coord$CoordList,
                                 WhichGrp_adjP.Separate=WhichGrp_adjP.Separate,
                                 FDRCutoff=obj$FDRCutoff,
                                 which_group_to_be_one=which_group_to_be_one,
                                 uniGroup=uniGroup)

  coord_0<-sep_checkSig$P_OR_group0
  coord_1<-sep_checkSig$P_OR_group1

  #################
  temp_ElogOR<-obj$P_OR_group0$ElogOR
  colN_temp<-colnames(temp_ElogOR)
  rowN_temp<-rownames(temp_ElogOR)


  Col_index<-coord_0$X
  Row_index<-nrow(temp_ElogOR)+1-coord_0$Y

  out<-data.frame(Gene1=rowN_temp[Row_index],
                  Gene2=colN_temp[Col_index],
                  coord_0[,c(3:6,9)],coord_1[,c(3:6,9)])
  colnames(out)[3:12]<-c(sprintf("%s.%s",c("LogOR","Raw.P","FDR.Within.Grp","FDR.All","Sig.In"),name_grp0),
                         sprintf("%s.%s",c("LogOR","Raw.P","FDR.Within.Grp","FDR.All","Sig.In"),name_grp1))
  row.names(out)<-1:nrow(out)

  ###########################
  if(is.null(WhichGrp_adjP.Separate)){
    out<-out[,c(1:4,6:9,11:12)]
  }else if(WhichGrp_adjP.Separate=="All"){
    out<-out[,c(1:4,6:9,11:12)]
  }else if(WhichGrp_adjP.Separate=="None"){
    out<-out[,c(1:4,7:9,12)]
  }else if (WhichGrp_adjP.Separate=="Each"){
    out<-out[,c(1:5,7:10,12)]
  }else if(WhichGrp_adjP.Separate==name_grp1){
    out<-out[,c(1:4,7:10,12)]
  }else if(WhichGrp_adjP.Separate==name_grp0){
    out<-out[,c(1:5,7:9,12)]
  }else{
      stop("The name to specify WhichGrp_adjP.Separate is not correct!")
    }

  return(out)
}




adjust_Single_P<-function(inputList=P_OR.list.single,
                          FDRCutoff=FDRCutoff,
                          adjust.method=adjust.method){
  coord_sepGroup<-extract_coordinate_singleGrp(dataListIn=inputList,mode="Single") #####
  combine_rawP<-c(coord_sepGroup[[1]]$rawP)
  combine_ElogOR<-c(coord_sepGroup[[1]]$ElogOR)
  combine_df<-data.frame(com_OR=combine_ElogOR,com_rawP=combine_rawP)
  combine_df$FDR<-p.adjust(combine_df$com_rawP,method=adjust.method)
  combine_df$scaledOR<-OR_Scale(OR_In=combine_df$com_OR)*0.97
  nRow_df<-nrow(combine_df)/2

  coord_sepGroup[[1]]$FDR_all<-combine_df$FDR
  coord_sepGroup[[1]]$OR_scaled<-combine_df$scaledOR

  inputList$CoordList<-coord_sepGroup
  inputList$combine_df<-combine_df
  inputList$FDRCutoff<-FDRCutoff
  return(inputList)
}




make_P_OR_Single<-function(mat=mutation_data,adjust.method=adjust.method){
  P_OR_single<-calculate_OD_P_oneGroup(dataIn=mat,adjust.method=adjust.method)
  return(P_OR_single)
}



adjust_SepGroup_P<-function(inputList=P_OR.list.sep.group,
                            FDRCutoff=FDRCutoff,
                            adjust.method=adjust.method){
  coord_sepGroup<-extract_coordinate_singleGrp(dataListIn=inputList,mode="Separate")
  combine_rawP<-c(coord_sepGroup[[1]]$rawP,coord_sepGroup[[2]]$rawP)
  combine_ElogOR<-c(coord_sepGroup[[1]]$ElogOR,coord_sepGroup[[2]]$ElogOR)
  combine_df<-data.frame(com_OR=combine_ElogOR,com_rawP=combine_rawP)
  combine_df$FDR<-p.adjust(combine_df$com_rawP,method=adjust.method)
  combine_df$scaledOR<-OR_Scale(OR_In=combine_df$com_OR)*0.97
  nRow_df<-nrow(combine_df)/2

  coord_sepGroup[[1]]$FDR_all<-combine_df$FDR[1:nRow_df]
  coord_sepGroup[[1]]$OR_scaled<-combine_df$scaledOR[1:nRow_df]
  coord_sepGroup[[2]]$FDR_all<-combine_df$FDR[(nRow_df+1):(2*nRow_df)]
  coord_sepGroup[[2]]$OR_scaled<-combine_df$scaledOR[(nRow_df+1):(2*nRow_df)]

  inputList$CoordList<-coord_sepGroup
  inputList$combine_df<-combine_df
  inputList$FDRCutoff<-FDRCutoff
  return(inputList)
}



make_P_OR_SeqGroup<-function(mat=mutation_data,
                             group=group,
                             which_group_to_be_one=which_group_to_be_one,
                             adjust.method=adjust.method){
  P_OR_group0<-calculate_OD_P_oneGroup(dataIn=mat[!group %in% which_group_to_be_one,],adjust.method=adjust.method)
  P_OR_group1<-calculate_OD_P_oneGroup(dataIn=mat[group %in% which_group_to_be_one,],adjust.method=adjust.method)
  out<-list(P_OR_group0=P_OR_group0,P_OR_group1=P_OR_group1)
  return(out)
}



make_table<-function(dataIn=tem){
  dataIn<-dataIn[!(dataIn[,1]=="" | dataIn[,1]=="NA" | is.na(dataIn[,1]) | dataIn[,2]=="" | dataIn[,2]=="NA" | is.na(dataIn[,2])),]

  F0S0<-sum(dataIn[,1]==0 & dataIn[,2]==0)
  F0S1<-sum(dataIn[,1]==0 & dataIn[,2]==1)
  F1S0<-sum(dataIn[,1]==1 & dataIn[,2]==0)
  F1S1<-sum(dataIn[,1]==1 & dataIn[,2]==1)
  tab_out<-matrix(c(F0S0,F0S1,F1S0,F1S1),nrow=2,byrow=T)
  row.names(tab_out)<-c("FirstCol_Zero","FirstCol_One")
  colnames(tab_out)<-c("SecondCol_Zero","SecondCol_One")
  return(tab_out)
}



calculate_OD_P_oneGroup<-function(dataIn=dataIn,adjust.method=adjust.method){
  n_var<-ncol(dataIn)
  ElogOR<-ElogOR_lower<-ElogOR_upper<-pvalue_SingleGroup<-matrix(Inf,nrow=n_var-1,ncol=n_var-1)
  row.names(ElogOR)<-row.names(ElogOR_lower)<-row.names(ElogOR_upper)<-row.names(pvalue_SingleGroup)<-colnames(dataIn)[1:(n_var-1)]
  colnames(ElogOR)<-colnames(ElogOR_lower)<-colnames(ElogOR_upper)<-colnames(pvalue_SingleGroup)<-colnames(dataIn)[2:n_var]

  for(i_first in 1:(n_var-1)){
    for(i_second in (i_first+1):n_var){
      tem<-dataIn[,c(i_first,i_second)]
      tab_tem<-make_table(dataIn=tem)
      fisherP_tem<-fisher.test(tab_tem)$p.value

      OD_tem=0;
      if(sum(tab_tem==0)>=1){tab_tem<-tab_tem+0.5;OD_tem<-tab_tem[1,1]*tab_tem[2,2]/(tab_tem[1,2]*tab_tem[2,1])}
      if(sum(tab_tem==0)==0){OD_tem<-tab_tem[1,1]*tab_tem[2,2]/(tab_tem[1,2]*tab_tem[2,1])}

      OD_log<-log(OD_tem);se_tem<-sqrt(1/tab_tem[1,1]+1/tab_tem[1,2]+1/tab_tem[2,1]+1/tab_tem[2,2])
      OD_low<-OD_log-qnorm(1-0.025)*se_tem
      OD_up<-OD_log+qnorm(1-0.025)*se_tem

      ElogOR[i_first,i_second-1]<-OD_log
      ElogOR_lower[i_first,i_second-1]<-OD_low
      ElogOR_upper[i_first,i_second-1]<-OD_up
      pvalue_SingleGroup[i_first,i_second-1]<-fisherP_tem
    }
  }
  which_notInf<-which(is.finite(pvalue_SingleGroup), arr.ind = TRUE)
  FDR_tem<-p.adjust(pvalue_SingleGroup[which_notInf],method=adjust.method)
  FDR_out<-pvalue_SingleGroup
  FDR_out[which_notInf]<-FDR_tem
  out<-list(ElogOR=ElogOR,
            ElogOR_lower=ElogOR_lower,
            ElogOR_upper=ElogOR_upper,
            rawP=pvalue_SingleGroup,
            P.adjusted.one.group=FDR_out)
  return(out)
}



coord_lowerTriangle<-function(size=1/240,cordx,cordy){
  if (size>1/4 | size<0) {
    stop("Size for coord_lowerTriangle is wrong (value should be in the range of 0~1/4)!")
  } #this value should be less than 1/4; small value will be in large plot
  poly1<-data.frame(x=c(size,1-size*3,size)+cordx-0.5, y=c(size,size,1-size*3)+cordy-0.5)
  return(poly1)
}



coord_upperTriangle<-function(size=1/240,cordx,cordy){
  if (size>1/4 | size<0) {
    stop("Size for coord_upperTriangle is wrong (value should be in the range of 0~1/4)!")
  }
  poly2<-data.frame(x=c(1-size,size*3,1-size)+cordx-0.5, y=c(1-size,1-size,size*3)+cordy-0.5)
  return(poly2)
}



legend_tri<-function(cordX=0,cordY=0,pos="upper"){
  size=0.125
  if(pos=="upper"){
    poly1<-data.frame(x=c(-0.5+size,0.5-size,0.5-size)+cordX, y=c(0.5-size,0.5-size,-0.5+size)+cordY)
  }else{
    poly1<-data.frame(x=c(-0.5+size,-0.5+size,0.5-size)+cordX, y=c(-0.5+size,0.5-size,-0.5+size)+cordY)
    }
  return(poly1)
}



check_NAInteraction<-function(dfIn=dataFrame){
  tab_df<-table(dfIn)
  NAInteraction<-"No"

  for(i in 1:2){
    tabIn<-tab_df[,,i]
    tab_ReLabel<-as.data.frame.matrix(tabIn)
    colnames(tab_ReLabel)<-c("First","Second")
    colum_zero<-apply(tab_ReLabel,2,function(x) sum(x==0))
    if(sum(colum_zero==2)>0){
      NAInteraction<-"Yes"
    }
  }

  tab_diffG<-table(dfIn[,c(2:3,1)])
  for(i in 1:2){
    tabIn<-tab_diffG[,,i]
    tab_ReLabel<-as.data.frame.matrix(tabIn)
    if(sum(tab_ReLabel)==0){
      NAInteraction<-"Yes"
    }
  }

  return(NAInteraction)
}


check_dropGrp<-function(dfIn=dataFrame){
  tab_df<-table(dfIn[,c(2:3,1)])
  dropGrp<-"No"

  for(i in 1:2){
    tabIn<-tab_df[,,i]
    tab_ReLabel<-as.data.frame.matrix(tabIn)
    if(sum(tab_ReLabel)==0){
      dropGrp<-"Yes"
    }
  }

  return(dropGrp)
}


assign_color<-function(OR_scaled=OR_scaled,
                       col_negOR=c("blue","light blue"),
                       col_zero="white",
                       col_posOR=c("#FDDBC7", "red")){
  len_OR_scaled<-length(OR_scaled)
  len_0<-sum(OR_scaled==0)
  len_p<-sum(OR_scaled>0)
  len_n<-sum(OR_scaled<0)
  #  len_NA<-sum(is.na(OR_In))

  if(len_0>0){ ####OR has zero
    if(len_p==0){####all Negative
      colsch<-c(colorRampPalette(c(col_negOR,c(rep(col_zero,len_0)))) (len_OR_scaled))
    }else if(len_n==0){####all Positive
      colsch<-c(colorRampPalette(c(rep(col_zero,len_0),len_p)) (len_OR_scaled))
    }else{####has positive and  negative
      colsch_neg0<-c(colorRampPalette(c(col_negOR,col_zero)) (len_n+1))
      colsch_neg<-colsch_neg0[-length(colsch_neg0)]

      colsch_pos0<-c(colorRampPalette(c(col_zero,col_posOR)) (len_p+1))
      colsch_pos<-colsch_pos0[-1]

      colsch<-c(colsch_neg,rep(col_zero,len_0),colsch_pos)
      #########
    }
  }else{ ####OR has no zero
    if(len_p==0){####all Negative
      colsch0<-c(colorRampPalette(c(col_negOR,col_zero)) (len_OR_scaled+1))
      colsch<-colsch0[-length(colsch0)]
    }else if(len_n==0){####all Positive
      colsch0<-c(colorRampPalette(c(col_zero,col_posOR)) (len_OR_scaled+1))
      colsch<-colsch0[-1]
    }else{####has positive and  negative
      colsch_neg0<-c(colorRampPalette(c(col_negOR,col_zero)) (len_n+1))
      colsch_neg<-colsch_neg0[-length(colsch_neg0)]

      colsch_pos0<-c(colorRampPalette(c(col_zero,col_posOR)) (len_p+1))
      colsch_pos<-colsch_pos0[-1]

      colsch<-c(colsch_neg,colsch_pos)
    }
  }
  col_plot<-colsch[match(OR_scaled,OR_scaled[order(OR_scaled)])]
  colList<-list(col_plot=col_plot,col_legend=colsch)
  return(colList)
}



OR_Scale<-function(OR_In=OR_combine){
  len_Elog<-length(OR_In)
  len_0<-sum(OR_In==0)
  len_p<-sum(OR_In>0)
  len_n<-sum(OR_In<0)
  #  len_NA<-sum(is.na(OR_In))

  ###scale data
  minE<-range(OR_In)[1]
  maxE<-range(OR_In)[2]
  datScal<-OR_In
  if(minE*maxE<0){
    datScal[which(OR_In<0)]<-OR_In[which(OR_In<0)]/abs(minE)
    datScal[which(OR_In>0)]<-OR_In[which(OR_In>0)]/abs(maxE)
  }else if(minE*maxE==0 & minE==0){
    datScal<-OR_In/abs(maxE)
  }else if(minE*maxE==0 & maxE==0){
    datScal<-OR_In/abs(minE)
  }else if(minE*maxE>0 & minE>0){
    datScal<-OR_In/abs(maxE)
  }else{
    datScal<-OR_In/abs(minE)
  }
  return(datScal)
}


extract_coordinate<-function(dataIn=dataIn,name="test"){
  finite_index<-which(is.finite(dataIn) | is.na(dataIn), arr.ind = TRUE)
  ##### transform to plot coordinate (nrow-row+1->y; col->X)
  coord_plot <- as.data.frame(finite_index)
  coord_plot[,1] <-  finite_index[,2]
  coord_plot[,2] <-  nrow(dataIn)+1-finite_index[,1]
  colnames(coord_plot)<-c("X","Y")
  coord_plot$value<-dataIn[finite_index]
  colnames(coord_plot)[3]<-name
  #coord_plot$FDR<-temp_FDR[finite_index]
  return(list(coord_plot))
}



###extract coordinate---single group
extract_coordinate_singleGrp<-function(dataListIn=inputList,mode=mode){
  if(mode=="Separate"){
    coord_SepGroup<-lapply(dataListIn,function(x){
      temp_ElogOR<-x$ElogOR
      finite_index<-which(is.finite(temp_ElogOR), arr.ind = TRUE)
      ##### transform to plot coordinate (nrow-row+1->y; col->X)
      coord_plot <- as.data.frame(finite_index)
      coord_plot[,1] <-  finite_index[,2]
      coord_plot[,2] <-  nrow(temp_ElogOR)+1-finite_index[,1]
      colnames(coord_plot)<-c("X","Y")
      coord_plot$ElogOR<-temp_ElogOR[finite_index]
      coord_plot$rawP<-x$rawP[finite_index]
      coord_plot$adjP_singleGroup<-x$P.adjusted.one.group[finite_index]
      return(coord_plot)
    })
  }
  if(mode=="Single"){
      temp_ElogOR<-dataListIn$ElogOR
      finite_index<-which(is.finite(temp_ElogOR), arr.ind = TRUE)
      ##### transform to plot coordinate (nrow-row+1->y; col->X)
      coord_plot <- as.data.frame(finite_index)
      coord_plot[,1] <-  finite_index[,2]
      coord_plot[,2] <-  nrow(temp_ElogOR)+1-finite_index[,1]
      colnames(coord_plot)<-c("X","Y")
      coord_plot$ElogOR<-temp_ElogOR[finite_index]
      coord_plot$rawP<-dataListIn$rawP[finite_index]
      coord_plot$adjP_singleGroup<-dataListIn$P.adjusted.one.group[finite_index]
      coord_SepGroup=list(Single=coord_plot)
  }
  return(coord_SepGroup)
}


####calculate ElogLOD and FDR#####
calculate_sig_interaction<-function(mat=gbmDat_F[,2:ncol(gbmDat_F)],
                                    group=gbmDat_F[,1],
                                    which_group_to_be_one="GBM_L",
                                    adjust.method="BH"){
#  library(logistf)
  bin_y<-ifelse(group %in% which_group_to_be_one,1,0)
  n_var<-ncol(mat)
  pvalue_DE<-OR_grp0<-OR_grp1<-pvalue_CO<-OR_CO<-matrix(Inf,nrow=n_var-1,ncol=n_var-1)
  row.names(pvalue_DE)<-row.names(OR_grp0)<-row.names(OR_grp1)<-row.names(pvalue_CO)<-row.names(OR_CO)<-colnames(mat)[1:(n_var-1)]
  colnames(pvalue_DE)<-colnames(OR_grp0)<-colnames(OR_grp1)<-colnames(pvalue_CO)<-colnames(OR_CO)<-colnames(mat)[2:n_var]

  for(i_first in 1:(n_var-1)){
    for(i_second in (i_first+1):n_var){
      print(sprintf("Evaluating DC P, pair1:%s; pair2:%s",i_first,i_second)) #evaluation interaction p

      df<-data.frame(grp=bin_y,mat[,c(i_first,i_second)])
      colnames(df)[c(2,3)]<-c("var1","var2")
      if(check_dropGrp(dfIn=df)=="Yes"){
        logistf_main<-logistf::logistf(var1~var2,data=df,firth=TRUE)
        main_p<-logistf_main$prob[c("var2")]
        main_expOR<-coef(logistf_main)[c("var2")]
        pvalue_CO[i_first,i_second-1]<-main_p
        OR_CO[i_first,i_second-1]<-main_expOR
      }else{
        logistf_main<-logistf::logistf(var1~grp+var2,data=df,firth=TRUE)
        main_p<-logistf_main$prob[c("var2")]
        main_expOR<-coef(logistf_main)[c("var2")]
        pvalue_CO[i_first,i_second-1]<-main_p
        OR_CO[i_first,i_second-1]<-main_expOR
      }


      if(check_NAInteraction(df)=="Yes"){
        pvalue_DE[i_first,i_second-1]<-OR_grp0[i_first,i_second-1]<-OR_grp1[i_first,i_second-1]<-NA  # not defined because of singularities
      }else{
        logistf_reg<-logistf::logistf(var1~grp*var2,data=df,firth=TRUE)
        interaction_p<-logistf_reg$prob[c("grp:var2")]
        pvalue_DE[i_first,i_second-1]<-interaction_p
        OR_grp0[i_first,i_second-1]<-coef(logistf_reg)["var2"]
        OR_grp1[i_first,i_second-1]<-coef(logistf_reg)["var2"]+coef(logistf_reg)["grp:var2"]
      }
    }
  }
  out<-list(pvalue_DE=pvalue_DE,
            OR_grp0=OR_grp0,
            OR_grp1=OR_grp1,
            pvalue_CO=pvalue_CO,
            OR_CO=OR_CO)
  return(out)
}



#########adjuste p value#####
adjust_inter_P<-function(inputList=interObj,
                         FDRCutoff=FDRCutoff,
                         adjust.method=adjust.method){
  pvalue_DE<-inputList$pvalue_DE
  pvalue_CO<-inputList$pvalue_CO

  ###adjust interaction p value first###########
  which_notInf<-which(is.finite(pvalue_DE), arr.ind = TRUE)
  FDR_tem<-p.adjust(pvalue_DE[which_notInf],method=adjust.method)
  FDR_inter_out<-pvalue_DE
  FDR_inter_out[which_notInf]<-FDR_tem

  sig_inter_coord<-which(!is.na(FDR_inter_out) & (FDR_inter_out<FDRCutoff), arr.ind = TRUE)
  ########main effect p value adjustement########
  which_main_notInf<-which(is.finite(pvalue_CO), arr.ind = TRUE)
  if(!is.na((sig_inter_coord)[1])){
    str_main<-sprintf("%s_%s",which_main_notInf[,1],which_main_notInf[,2])
    str_inter<-sprintf("%s_%s",sig_inter_coord[,1],sig_inter_coord[,2])

    which_main_notInf_notInterSig<-which_main_notInf[!str_main %in% str_inter,]
    FDR_main_adj<-p.adjust(pvalue_CO[which_main_notInf_notInterSig],method=adjust.method)

    FDR_CO<-pvalue_CO
    FDR_CO[sig_inter_coord]<-NA  ## set this value to 1
    FDR_CO[which_main_notInf_notInterSig]<-FDR_main_adj
  }else{
    FDR_main_adj<-p.adjust(pvalue_CO[which_main_notInf],method=adjust.method)
    FDR_CO<-pvalue_CO
    FDR_CO[which_main_notInf]<-FDR_main_adj
  }
  inputList$FDR_DE=FDR_inter_out
  inputList$FDR_CO=FDR_CO
  inputList$FDRCutoff=FDRCutoff
  return(inputList)
}




get_coordList<-function(obj=obj,n=n){
  CoordList<-mapply(extract_coordinate,obj[1:n],names(obj)[1:n])
  return(CoordList)
}




#########df_coord_P_OD#####
make_Coord_df<-function(obj=list.P.adj,
                        DE_mode=DE_mode,
                        FDRCutoff=0.05,
                        col_negOR=c("blue","light blue"),
                        col_zero="white",
                        col_posOR=c("pink", "red")){
  if(DE_mode=="DC"){
    CoordList<-get_coordList(obj=obj,n=7)
    sig_inter_index<-which(CoordList$FDR_DE$FDR_D<FDRCutoff)
    length_sigInter<-length(sig_inter_index)

    SigOR0<-CoordList$OR_grp0[sig_inter_index,"OR_grp0"]
    SigOR1<-CoordList$OR_grp1[sig_inter_index,"OR_grp1"]

    #####Main effect#####
    seq_nrow<-1:nrow(CoordList$OR_CO)
    sig_CO_index<-seq_nrow[!seq_nrow %in% sig_inter_index]
    length_main<-length(sig_CO_index)
    OR_Main<-CoordList$OR_CO[sig_CO_index,"OR_CO"]

    ##combine OR#######
    OR_combine<-c(SigOR0,SigOR1,OR_Main);length(OR_combine)
    OR_scaled<-OR_Scale(OR_In=OR_combine)*0.97

    ###assign scale data into list###
    CoordList$OR_grp0$ScaledOR<-CoordList$OR_grp1$ScaledOR<-CoordList$OR_CO$ScaledOR<-NA
    CoordList$OR_grp0$ScaledOR[sig_inter_index]<-OR_scaled[1:length_sigInter]
    CoordList$OR_grp1$ScaledOR[sig_inter_index]<-OR_scaled[(length_sigInter+1):(2*length_sigInter)]
    CoordList$OR_CO$ScaledOR[sig_CO_index]<-OR_scaled[(2*length_sigInter+1):length(OR_scaled)]
    ####index and scaled color###
    cols<-assign_color(OR_scaled=OR_scaled,col_negOR=col_negOR,col_zero=col_zero,col_posOR=col_posOR)

    CoordList$OR_grp0$Color<-CoordList$OR_grp1$Color<-CoordList$OR_CO$Color<-NA
    CoordList$OR_grp0$Color[sig_inter_index]<-cols$col_plot[1:length_sigInter]
    CoordList$OR_grp1$Color[sig_inter_index]<-cols$col_plot[(length_sigInter+1):(2*length_sigInter)]
    CoordList$OR_CO$Color[sig_CO_index]<-cols$col_plot[(2*length_sigInter+1):length(cols$col_plot)]

    OR_all<-list(OR_combine=OR_combine,OR_scaled=OR_scaled)
    out<-list(matrix_obj=obj,CoordList=CoordList,colors=cols,OR_all=OR_all)
  }
  if(DE_mode=="Separate"){
    CoordList<-obj$CoordList
    OR_scaled<-obj$combine_df$scaledOR
    cols<-assign_color(OR_scaled=OR_scaled,col_negOR=col_negOR,col_zero=col_zero,col_posOR=col_posOR)
    nrow_data<-nrow(CoordList$P_OR_group0)
    CoordList$P_OR_group0$Color<-cols$col_plot[1:nrow_data]
    CoordList$P_OR_group1$Color<-cols$col_plot[(nrow_data+1):(2*nrow_data)]

    OR_all<-list(OR_combine=obj$combine_df$com_OR,OR_scaled=OR_scaled)
    out<-list(matrix_obj=obj,CoordList=CoordList,colors=cols,OR_all=OR_all)
  }
  if(DE_mode=="Single"){
    CoordList<-obj$CoordList
    OR_scaled<-obj$combine_df$scaledOR
    cols<-assign_color(OR_scaled=OR_scaled,col_negOR=col_negOR,col_zero=col_zero,col_posOR=col_posOR)
    CoordList$Single$Color<-cols$col_plot

    OR_all<-list(OR_combine=obj$combine_df$com_OR,OR_scaled=OR_scaled)
    out<-list(matrix_obj=obj,CoordList=CoordList,colors=cols,OR_all=OR_all)
  }
  return(out)
}




xylim<-function(rang_X=rang_X,
                rang_Y=rang_Y,
                sig_CO.FDRCutoff=sig_CO.FDRCutoff,
                max_length_geneName=max_length_geneName,
                label.gene.offset.vertical=label.gene.offset.vertical,
                label.gene.offset.horiz=label.gene.offset.horiz,
                uniGroup=uniGroup,
                label.gene.cex=label.gene.cex,
                legend.offset.X=legend.offset.X,
                legend.width=legend.width,
                legend.text.cex=legend.text.cex){
  unigroupLength<-max(strwidth(uniGroup,cex=label.gene.cex,units="inches"))
  FDRLength<-max(strwidth(sprintf("FDR_DC<%s",sig_CO.FDRCutoff),cex=label.gene.cex,units="inches"),
                 strwidth(sprintf("FDR_CO<%s",sig_CO.FDRCutoff),cex=label.gene.cex,units="inches"))
  colLegendWidth<-strwidth(as.character(-99.99),cex=legend.text.cex,units="inches")+legend.width+legend.offset.X+legend.offset.X+0.2+1



  max_X_right_Width<-max(unigroupLength+1.5+label.gene.offset.horiz+1, ####labels of group
                         FDRLength+1.5+label.gene.offset.horiz+1,      #### labels of FDR
                         colLegendWidth)                                 ##### labels of col_legend
  max_X_left_Width<-label.gene.offset.horiz+max_length_geneName+1
  max_y_right_Width<-label.gene.offset.vertical+max_length_geneName+1

  ylim_left<-rang_Y[1]-1
  ylim_right<-rang_Y[2]+max_y_right_Width
  xlim_left<-rang_X[1]-max_X_left_Width
  xlim_right<-rang_X[2]+max_X_right_Width
  out<-list(xlim_left=xlim_left,
            xlim_right=xlim_right,
            ylim_left=ylim_left,
            ylim_right=ylim_right)
  return(out)
}




check_Sig_SepGrp<-function(dataIn=obj_coord$CoordList,
                           WhichGrp_adjP.Separate="All",
                           FDRCutoff=FDRCutoff,
                           which_group_to_be_one=which_group_to_be_one,
                           uniGroup=uniGroup){
  group_0_name<-uniGroup[!uniGroup %in% which_group_to_be_one]

  if(is.null(WhichGrp_adjP.Separate)){
    list_signSig<-lapply(dataIn,function(x){x$Sig<-ifelse(x$FDR_all<FDRCutoff,"Yes","No");return(x)})
  }else if(WhichGrp_adjP.Separate=="All"){
    list_signSig<-lapply(dataIn,function(x){x$Sig<-ifelse(x$FDR_all<FDRCutoff,"Yes","No");return(x)})
  }else if(WhichGrp_adjP.Separate=="None"){
    list_signSig<-lapply(dataIn,function(x){x$Sig<-ifelse(x$rawP<FDRCutoff,"Yes","No");return(x)})
  }else if(WhichGrp_adjP.Separate=="Each"){
    list_signSig<-lapply(dataIn,function(x){x$Sig<-ifelse(x$adjP_singleGroup<FDRCutoff,"Yes","No");return(x)})
  }else if(WhichGrp_adjP.Separate==which_group_to_be_one){
    list_signSig<-dataIn
    list_signSig$P_OR_group1$Sig<-ifelse(list_signSig$P_OR_group1$adjP_singleGroup<FDRCutoff,"Yes","No")
    list_signSig$P_OR_group0$Sig<-ifelse(list_signSig$P_OR_group0$rawP<FDRCutoff,"Yes","No")
  }else if(WhichGrp_adjP.Separate==group_0_name){
    list_signSig<-dataIn
    list_signSig$P_OR_group1$Sig<-ifelse(list_signSig$P_OR_group1$rawP<FDRCutoff,"Yes","No")
    list_signSig$P_OR_group0$Sig<-ifelse(list_signSig$P_OR_group0$adjP_singleGroup<FDRCutoff,"Yes","No")
  }else{
    stop("The name to specify WhichGrp_adjP.Separate is not correct!")
  }
  return(list_signSig)
}





############# add gene name to the plot########
add_geneName<-function(coordIn_main=coordIn_main,
                       label.gene.offset.vertical=label.gene.offset.vertical,
                       label.gene.offset.horiz=label.gene.offset.horiz,
                       plot_xlab=plot_xlab,
                       plot_ylab=plot_ylab,
                       label.gene.col=label.gene.col,
                       label.gene.font=label.gene.font,
                       label.gene.cex=label.gene.cex){
  pos.xlabel <- cbind(min(coordIn_main[,1]):max(coordIn_main[,1]), max(coordIn_main[,2])+label.gene.offset.vertical)
  pos.ylabel <- cbind(min(coordIn_main[,1]):max(coordIn_main[,1])-label.gene.offset.horiz, max(coordIn_main[,2]):min(coordIn_main[,2]))

  text(pos.xlabel[,1], pos.xlabel[,2], plot_xlab, srt=90,adj=0,col=label.gene.col, font=label.gene.font, cex =label.gene.cex,xpd=T)
  text(pos.ylabel[,1], pos.ylabel[,2], plot_ylab,col =label.gene.col,font=label.gene.font,cex=label.gene.cex,pos=2, adj=0,xpd=T)
}




###########add symbols to the plot########
add_symbols<-function(DE_mode=DE_mode,
                      coordIn_main=coordIn_main,
                      coordIn_FDRmain_sig=coordIn_FDRmain_sig,
                      sig_CO.symbol=sig_CO.symbol,
                      sig_CO.cex=sig_CO.cex,
                      sig_CO.col=sig_CO.col,
                      coordIn_OR_grp0=coordIn_OR_grp0,
                      coordIn_OR_grp1=coordIn_OR_grp1,
                      sig_DE.col=sig_DE.col,
                      sig_DE.lwd=sig_DE.lwd,
                      triangle.lwd=triangle.lwd,
                      grid.color=grid.color,
                      grid.lwd=grid.lwd,
                      FDRCutoff=FDRCutoff){
  if(DE_mode=="DC"){
    symbols(coordIn_main[,1:2], add = TRUE, inches = FALSE,
            rectangles = matrix(rep(sqrt(abs(coordIn_main$ScaledOR)),each=2), nrow = nrow(coordIn_main), ncol = 2,byrow=T),
            bg =coordIn_main$Color, fg =coordIn_main$Color,lwd=0.1)
    text(coordIn_FDRmain_sig[,1],coordIn_FDRmain_sig[,2],labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col)
    ##sig iner
    if(!is.na(coordIn_OR_grp0[1,1])){
      symbols(coordIn_OR_grp0[,1:2], add = TRUE, inches = FALSE,rectangles = matrix(1, nrow = nrow(coordIn_OR_grp0), ncol = 2),
              bg =NA, fg=sig_DE.col,lwd=sig_DE.lwd)

      coord_upper<-coordIn_OR_grp1;coord_lower<-coordIn_OR_grp0
      for(i in 1:nrow(coordIn_OR_grp0)){
        segments(coordIn_OR_grp0[i,1]-0.5,coordIn_OR_grp0[i,2]+0.5,coordIn_OR_grp0[i,1]+0.5,coordIn_OR_grp0[i,2]-0.5,col=sig_DE.col,lwd=sig_DE.lwd/2)
        upperCoord<-coord_upperTriangle(size=1/4-sqrt(abs(coord_upper$ScaledOR[i]))/4,cordx=coord_upper$X[i],cordy=coord_upper$Y[i])
        lowerCoord<-coord_lowerTriangle(size=1/4-sqrt(abs(coord_lower$ScaledOR[i]))/4,cordx=coord_lower$X[i],cordy=coord_lower$Y[i])
        polygon(upperCoord$x, upperCoord$y, xpd = NA, col = coord_upper$Color[i],lty = 1, lwd = triangle.lwd, border = NA)
        polygon(lowerCoord$x, lowerCoord$y, xpd = NA, col = coord_lower$Color[i],lty = 1, lwd = triangle.lwd, border = NA)
      }
    }
  }

  if(DE_mode=="Separate"){
    coord_upper<-coordIn_OR_grp1;coord_lower<-coordIn_OR_grp0
    for(i in 1:nrow(coordIn_OR_grp0)){
      segments(coordIn_OR_grp0[i,1]-0.5,coordIn_OR_grp0[i,2]+0.5,coordIn_OR_grp0[i,1]+0.5,coordIn_OR_grp0[i,2]-0.5,col=grid.color,lwd=grid.lwd/10)
      upperCoord<-coord_upperTriangle(size=1/4-sqrt(abs(coord_upper$OR_scaled[i]))/4,cordx=coord_upper$X[i],cordy=coord_upper$Y[i])
      lowerCoord<-coord_lowerTriangle(size=1/4-sqrt(abs(coord_lower$OR_scaled[i]))/4,cordx=coord_lower$X[i],cordy=coord_lower$Y[i])
      polygon(upperCoord$x, upperCoord$y, xpd = NA, col = coord_upper$Color[i],lty = 1, lwd = triangle.lwd, border = NA)
      polygon(lowerCoord$x, lowerCoord$y, xpd = NA, col = coord_lower$Color[i],lty = 1, lwd = triangle.lwd, border = NA)
    }
    if(sum(coord_upper$Sig=="Yes")>0){
      sig_upper_label<-coord_upper[coord_upper$Sig=="Yes",]
      text(sig_upper_label[,1]+0.25,sig_upper_label[,2]+0.25,labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col)}
    if(sum(coord_lower$Sig=="Yes")>0){
      sig_lower_label<-coord_lower[coord_lower$Sig=="Yes",]
      text(sig_lower_label[,1]-0.25,sig_lower_label[,2]-0.25,labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col)}
  }

  if(DE_mode=="Single"){
    symbols(coordIn_main[,1:2], add = TRUE, inches = FALSE,
            rectangles = matrix(rep(sqrt(abs(coordIn_main$OR_scaled)),each=2), nrow = nrow(coordIn_main), ncol = 2,byrow=T),
            bg =coordIn_main$Color, fg =coordIn_main$Color,lwd=0.1)
    coordIn_FDRmain_sig_Single<-coordIn_main[coordIn_main$FDR_all<FDRCutoff,]
    if(nrow(coordIn_FDRmain_sig_Single)>0){
      text(coordIn_FDRmain_sig_Single[,1],coordIn_FDRmain_sig_Single[,2],labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col)
    }
  }
}


#################################
add_legend<-function(coordIn_main=coordIn_main,
                     legend.offset.X=legend.offset.X,
                     legend.triangle.position=legend.triangle.position,
                     legend.y.below.manualset=legend.y.below.manualset,
                     legend.y.height.manualset=legend.y.height.manualset,
                     obj_coord=obj_coord,
                     legend.width=legend.width,
                     legend.Height=legend.Height,
                     text.logor.offset=text.logor.offset,
                     legend.text.cex=legend.text.cex,
                     legend.text.length=legend.text.length){
  legend.x=max(coordIn_main[,1]+legend.offset.X)
  if(is.null(legend.triangle.position)){legend.triangle.position<-ifelse(length(unique(coordIn_main[,2]))>10,"Right","LeftBottom")}

  if(is.null(legend.y.below.manualset)){
    diff.y.range<-diff(range(coordIn_main[,2]))
    if(legend.triangle.position=="Right"){
      legend.y<-quantile(unique(coordIn_main[,2]),0.4);legend.Height<-quantile(unique(coordIn_main[,2]),0.85)-legend.y
      if(diff.y.range<=6){
        legend.y<-quantile(unique(coordIn_main[,2]),0.6);legend.Height<-max(coordIn_main[,2])-legend.y+0.5
      }
    }else{
      legend.y<-quantile(unique(coordIn_main[,2]),0.25);legend.Height<-quantile(unique(coordIn_main[,2]),0.75)-legend.y
    }
  }else{
    legend.y<-legend.y.below.manualset;legend.Height<-legend.y.height.manualset
  }
  df_Legend<-data.frame(colLegend=obj_coord$colors$col_legend,OR=obj_coord$OR_all$OR_combine[order(obj_coord$OR_all$OR_combine)])

  plotrix::color.legend(xl=legend.x, yb=legend.y, xr=legend.x+legend.width, yt =legend.y+legend.Height, legend = "" , gradient="y", rect.col=df_Legend$colLegend, align="rb")
  if(is.null(text.logor.offset)){text.logor.offset=max(coordIn_main[,2])/40}
  text(legend.x,legend.y+legend.Height+text.logor.offset,"LogOR",adj=0,cex=legend.text.cex,xpd=T)
#

  ##legend text---coord
  legend.text.x<-legend.x+legend.width+0.2
  legend.text.pos.y<-seq(legend.y,legend.y+legend.Height,length=legend.text.length)
  legend.text.lab<-seq(range(df_Legend$OR)[1],range(df_Legend$OR)[2],length=legend.text.length)
  df.legend.text<-data.frame(x=legend.text.x,y=legend.text.pos.y,lab=as.character(round(legend.text.lab,2)))
  text(df.legend.text$x,df.legend.text$y,df.legend.text$lab,adj=0,cex=legend.text.cex,xpd=T)

  return(legend.y)
}


###########################
add_triangle_FDR_LeftBottom<-function(coordIn_main=coordIn_main,
                                      DE_mode=DE_mode,
                                      coordIn_OR_grp0=coordIn_OR_grp0,
                                      coordIn_OR_grp1=coordIn_OR_grp1,
                                      triangle.lwd=triangle.lwd,
                                      label.gene.offset.horiz=label.gene.offset.horiz,
                                      uniGroup=uniGroup,
                                      which_group_to_be_one=which_group_to_be_one,
                                      FDRCutoff=FDRCutoff,
                                      sig_DE.col=sig_DE.col,
                                      sig_DE.lwd=sig_DE.lwd,
                                      sig_CO.FDRCutoff=sig_CO.FDRCutoff,
                                      sig_CO.symbol=sig_CO.symbol,
                                      sig_CO.cex=sig_CO.cex,
                                      sig_CO.col=sig_CO.col,
                                      label.gene.cex=label.gene.cex,
                                      legend.text.cex=legend.text.cex,
                                      coordIn_FDRmain_sig=coordIn_FDRmain_sig){

  min.coordX<-min(coordIn_main[,1]);min.coordY<-min(coordIn_main[,2])
  ###############
  if(!is.na(coordIn_OR_grp0[1,1]) & !DE_mode=="Single"){#there is significant DE items or in sepearte mode
    legend_upperCoord<-legend_tri(cordX=min.coordX,cordY=min.coordY+1,pos="upper")
    legend_lowerCoord<-legend_tri(cordX=min.coordX,cordY=min.coordY,pos="lower")
    polygon(legend_upperCoord$x, legend_upperCoord$y, xpd = NA, col ="black",lty = 1, lwd = triangle.lwd, border = NA)
    polygon(legend_lowerCoord$x, legend_lowerCoord$y, xpd = NA, col = "black",lty = 1, lwd = triangle.lwd, border = NA)
    text(min.coordX+label.gene.offset.horiz,min.coordY,sprintf("%s",uniGroup[!uniGroup %in% which_group_to_be_one]),adj=0,cex=label.gene.cex,xpd=T)
    text(min.coordX+label.gene.offset.horiz,min.coordY+1,sprintf("%s",which_group_to_be_one),adj=0,cex=label.gene.cex,xpd=T)
    ##the above code adds polygon legend as well as the labeles.

    ##### for DE mode--add sig interaction labels; For sep.analysis, add FDR<0.05########
    if(DE_mode=="DC"){
      symbols(x=min.coordX,y=min.coordY+2, add = TRUE, inches = FALSE, rectangles = matrix(0.8, nrow =1, ncol = 2,byrow=T),bg =NA, fg=sig_DE.col,lwd=sig_DE.lwd)
      segments(min.coordX-0.4,min.coordY+2.4,min.coordX+0.4,min.coordY+1.6,col=sig_DE.col,lwd=sig_DE.lwd/2)
      text(min.coordX+label.gene.offset.horiz,min.coordY+2,sprintf("FDR_DC<%s",sig_CO.FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)
      #the above code adds significant square labels

      if(sum(!is.na(coordIn_FDRmain_sig[,1]))>0){ ###
        text(min.coordX,min.coordY+3,labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col,xpd=T)
        text(min.coordX+label.gene.offset.horiz,min.coordY+3,sprintf("FDR_CO<%s",sig_CO.FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)
      }#the above code adds FDR_CO<0.05
    }

    if(DE_mode=="Separate"){### end of "if" of DE_mode
      coord_upper<-coordIn_OR_grp1;coord_lower<-coordIn_OR_grp0
      if(sum(coord_lower$Sig=="Yes")>0 | sum(coord_upper$Sig=="Yes")>0){#for seperate analysis group
        text(min.coordX,min.coordY+2,labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col,xpd=T)
        text(min.coordX+label.gene.offset.horiz,min.coordY+2,sprintf("FDR<%s",FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)
      }#end of "for seperate analysis group"
    }## end of DE_mode=="Sep"
    #####"for DE mode--add sig interaction labels; For seperate group analysis, add FDR<0.05"########

  }else{###end of "if(!is.na(coordIn_OR_grp0[1,1]))";  check if main in DE mode has FDR<0.05
    if(DE_mode=="DC"){
      if(sum(!is.na(coordIn_FDRmain_sig[,1]))>0){ ###no sig interaction, but sig main effect
        text(min.coordX,min.coordY,labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col,xpd=T)
        text(min.coordX+label.gene.offset.horiz,min.coordY,sprintf("FDR<%s",sig_CO.FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)
      }### end of "no sig interaction, but sig main effect"
    }
    if(DE_mode=="Single"){
      coordIn_FDRmain_sig_Single<-coordIn_main[coordIn_main$FDR_all<FDRCutoff,]
      if(nrow(coordIn_FDRmain_sig_Single)>0){ ###Single mode
        text(min.coordX,min.coordY,labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col,xpd=T)
        text(min.coordX+label.gene.offset.horiz,min.coordY,sprintf("FDR<%s",FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)
      }###
    }
  }
}


################################
add_triangle_FDR_Right<-function(legend.y=legend.y.returned,
                                 coordIn_main=coordIn_main,
                                 DE_mode=DE_mode,
                                 triangle.lwd=triangle.lwd,
                                 label.gene.offset.horiz=label.gene.offset.horiz,
                                 uniGroup=uniGroup,
                                 which_group_to_be_one=which_group_to_be_one,
                                 legend.text.cex=legend.text.cex,
                                 legend.offset.X=legend.offset.X,
                                 legend.y.below.manualset=legend.y.below.manualset,
                                 sig_DE.col=sig_DE.col,
                                 sig_DE.lwd=sig_DE.lwd,
                                 sig_CO.FDRCutoff=sig_CO.FDRCutoff,
                                 label.gene.cex=label.gene.cex,
                                 coordIn_OR_grp0=coordIn_OR_grp0,
                                 coordIn_OR_grp1=coordIn_OR_grp1,
                                 FDRCutoff=FDRCutoff,
                                 sig_CO.symbol=sig_CO.symbol,
                                 sig_CO.cex=sig_CO.cex,
                                 sig_CO.col=sig_CO.col,
                                 coordIn_FDRmain_sig=coordIn_FDRmain_sig){
  min.coordX<-max(coordIn_main[,1])+legend.offset.X+0.5

  min.coordY_0<-legend.y-2;min.coordY_1<-quantile(unique(coordIn_main[,2]),0.25)
  min.coordY<-max(min.coordY_0,min.coordY_1)
  if(!is.null(legend.y.below.manualset)){min.coordY<-legend.y.below.manualset-2}


  rangeDiff_Y<-diff(range(coordIn_main[,2]))
  triangleSize=0.05;min.coordY.seq<-min.coordY-1*0:3;RecWH<-0.8
  if(rangeDiff_Y<12){triangleSize=0.125;min.coordY.seq<-min.coordY-0.6*0:3;RecWH<-0.5}
  diag.set.value<-RecWH/2

  if(!is.na(coordIn_OR_grp0[1,1]) & !DE_mode=="Single"){
    legend_upperCoord<-coord_upperTriangle(size=triangleSize,cordx=min.coordX,cordy=min.coordY.seq[1])
    legend_lowerCoord<-coord_lowerTriangle(size=triangleSize,cordx=min.coordX,cordy=min.coordY.seq[2])
    polygon(legend_upperCoord$x-0.125, legend_upperCoord$y-0.125, xpd = NA, col ="black",lty = 1, lwd = triangle.lwd, border = NA)
    polygon(legend_lowerCoord$x+0.125, legend_lowerCoord$y+0.125, xpd = NA, col = "black",lty = 1, lwd = triangle.lwd, border = NA)

    text(min.coordX+label.gene.offset.horiz,min.coordY.seq[1],sprintf("%s",which_group_to_be_one),adj=0,cex=legend.text.cex,xpd=T)
    text(min.coordX+label.gene.offset.horiz,min.coordY.seq[2],sprintf("%s",uniGroup[!uniGroup %in% which_group_to_be_one]),adj=0,cex=legend.text.cex,xpd=T)

    #####
    if(DE_mode=="DC"){
      #
      symbols(x=min.coordX,y=min.coordY.seq[3], add = TRUE, inches = FALSE,
              rectangles = matrix(RecWH, nrow =1, ncol = 2,byrow=T),bg =NA, fg=sig_DE.col,lwd=sig_DE.lwd)
      segments(min.coordX-diag.set.value,min.coordY.seq[3]+diag.set.value,min.coordX+diag.set.value,min.coordY.seq[3]-diag.set.value,col=sig_DE.col,lwd=sig_DE.lwd/2)
      text(min.coordX+label.gene.offset.horiz,min.coordY.seq[3],sprintf("FDR_DC<%s",sig_CO.FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)
      #
      if(sum(!is.na(coordIn_FDRmain_sig[,1]))>0){
        text(min.coordX,min.coordY.seq[4],labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col,xpd=T)
        text(min.coordX+label.gene.offset.horiz,min.coordY.seq[4],sprintf("FDR_CO<%s",sig_CO.FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)
      }
    }

    if(DE_mode=="Separate"){
      #####
      coord_upper<-coordIn_OR_grp1;coord_lower<-coordIn_OR_grp0
      if(sum(coord_lower$Sig=="Yes")>0 | sum(coord_upper$Sig=="Yes")>0){#
        text(min.coordX,min.coordY.seq[3],labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col,xpd=T)
        text(min.coordX+label.gene.offset.horiz,min.coordY.seq[3],sprintf("FDR<%s",FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)}#
    }



  }else{########
    if(DE_mode=="DC"){
      if(sum(!is.na(coordIn_FDRmain_sig[,1]))>0){ ###no sig interaction, but sig main effect
        text(min.coordX-0.5,min.coordY.seq[1],labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col,xpd=T)
        text(min.coordX+label.gene.offset.horiz-0.5,min.coordY.seq[1],sprintf("FDR<%s",sig_CO.FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)
      }### end of "no sig interaction, but sig main effect"
    }
    if(DE_mode=="Single"){
      coordIn_FDRmain_sig_Single<-coordIn_main[coordIn_main$FDR_all<FDRCutoff,]
      if(nrow(coordIn_FDRmain_sig_Single)>0){ ###Single mode
        text(min.coordX-0.5,min.coordY.seq[1],labels=sig_CO.symbol,cex=sig_CO.cex,col=sig_CO.col,xpd=T)
        text(min.coordX+label.gene.offset.horiz-0.5,min.coordY.seq[1],sprintf("FDR<%s",FDRCutoff),adj=0,cex=legend.text.cex,xpd=T)
      }###
    }

  }########
}##end of Right legend





