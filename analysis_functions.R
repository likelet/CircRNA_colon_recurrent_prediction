#data processing functions
suppressMessages(library("sva"))
suppressMessages(library("impute"))
suppressMessages(library("data.table"))
#lasso functions
suppressMessages(library("glmnet"))
#surival analysis
suppressMessages(library("survival"))
#plot functions
suppressMessages(library("ggplot2"))
suppressMessages(library("survminer"))
suppressMessages(library("hdnom"))
suppressMessages(library("papeR"))
suppressMessages(library("ggsci"))
suppressMessages(library("grid"))
suppressMessages(library("ggthemes"))
suppressMessages(library("cowplot"))
suppressMessages(library("papeR"))
suppressMessages(library("pheatmap"))
#perfunmance evaluation functions
suppressMessages(library("risksetROC"))
suppressMessages(library("pROC"))
suppressMessages(library("ROCR"))
suppressMessages(library("gridExtra"))
suppressMessages(library("rms"))
batchUnivarCOXfun<-function(surv,data){
  df<-data.frame(name=c("coef","Hazard Ratio","CI (lower)","CI (upper)","se(coef)","z","Pr(>|Z|)","C-index","C-index-se","testCph.p"))
  pro=ncol(data)/10000
  j=1;
  k=1;
  for (i in 1:ncol(data)) {
    fml=as.formula(paste("surv ~ ",names(data)[i]))
    fit<-coxph(surv ~ data[,i])
    p<-data.frame(summary(fit)$coefficients)[1,5]
    test=cox.zph(coxph(surv ~ data[,i]))
    testp= test$table[3]
    #     x<-c(x,y)
    #     coefficients=summary(a)$coefficients
    y=prettify(summary(fit))
    y=y[,c(-1,-8,-9)]
    x=c(t(y)[,1],p,summary(fit)$concordance,testp)
    df=data.frame(df,x)
    if(i==pro*k){
      print(paste(j*k,"items finished"))
      k=k+1
    }
  }
  row.names(df)<-df[,1]
  df<-t(df[,-1])
  row.names(df)=names(data)
  return(df)
}
batchTtest <- function(df,class,rownumber=2){
  
  singleTtest=function(value,class){
    classV=unique(class)
    x=value[class==classV[1]]
    y=value[class==classV[2]]
    
    return(wilcox.test(x,y)$p.value)
  }
  singlefoldchange=function(value,class){
    
    classV=unique(class)
    x=value[class==classV[1]]
    y=value[class==classV[2]]
    
    return(mean(x)-mean(y))
  }
  valuelist=apply(df,rownumber,singleTtest,class=class)
  fclist=apply(df,rownumber,singlefoldchange,class=class)
  fdrlist=p.adjust(valuelist)
  
  nm=c()
  if(rownumber==2){
    nm=colnames(df)
  }else{
    nm=row.names(df)
  }
  return(data.frame(ID=nm,fc=fclist,Pvalue=valuelist,fdr=fdrlist))
}
plotriskscoreInpatient=function(riskscorevec,status){
  df=data.frame(Riskscore=riskscorevec,Isalive=as.factor(status))
  df=df[with(df, order(-riskscorevec)), ]
  df=cbind(patientsnumber=seq(1:nrow(df)),df)
  df$patientsnumber=factor(df$patientsnumber,levels=seq(1:nrow(df)))
  p=ggplot(df,aes(x=patientsnumber,y=Riskscore,fill=as.factor(Isalive)))+geom_bar(stat="identity")+ggtitle(deparse(substitute(riskscorevec)))
  p=p+theme_classic()+theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            legend.position=c(.9,.8))+
    scale_fill_manual(values=c("#56B4E9", "#E69F00"), 
                      name="")
  return(p)
}
#cut vector by median
cutbymedian=function(vec){
  md=median(as.numeric(vec))
  vec2=rep("High",length(vec))
  vec2[vec<md]="Low"
  return(as.factor(vec2))
}
cutbymedian_df=function(df){
  df2=data.frame(row.names(df))
  for(i in 1:ncol(df)){
    df2= cbind(df2,cutbymedian(df[,i]))
  }
  df2=df2[,-1]
  colnames(df2)=colnames(df)
  return(df2)
  
}
cugbymedian=function(vec){
  tempx=rep(1,length(vec))
  tempx[vec<median(vec)]=0
  vec=tempx
  return(vec)
  
}
batchWilcox=function(df,class,rownumber=2){
  singleWilcox=function(value,class){
    classV=unique(class)
    x=value[class==classV[1]]
    y=value[class==classV[2]]
    
    return(wilcox.test(x,y)$p.value)
  }
  valuelist=apply(df,rownumber,singleWilcox,class=class)
  
  nm=c()
  if(rownumber==2){
    nm=colnames(df)
  }else{
    nm=row.names(df)
  }
  return(data.frame(ID=nm,Pvalue=valuelist))
}
batchfoldchange=function(df,class,rownumber=2){
  singleWilcox=function(value,class){
    classV=unique(class)
    x=value[class==classV[1]]
    y=value[class==classV[2]]
    
    return(mean(y)/mean(x))
  }
  valuelist=apply(df,rownumber,singleWilcox,class=class)
  print(paste(unique(class)[2],"/",unique(class)[1]))
  nm=c()
  if(rownumber==2){
    nm=colnames(df)
  }else{
    nm=row.names(df)
  }
  return(data.frame(ID=nm,foldchange=valuelist))
}

#plot function
savePlot <- function(myPlot) {
  pdf(paste(deparse(substitute(myPlot)),".pdf",sep=""))
  print(myPlot)
  dev.off()
}
#write out table ----
zhaoqi.write.csv.survival <- function(analysisDF){
  analysisDF$HRci=paste(analysisDF$`Hazard Ratio`," (",analysisDF$`CI (lower)`,"-",analysisDF$`CI (upper)`,")",sep="")
  df<-analysisDF[,c(1,2,10,7,8)]
  
  return(df)
}
#confusing matrix 
#preloading functions 
#from https://stackoverflow.com/questions/24716337/test-quality-of-logistic-regression-model-using-confusion-matrix-and-roc-curve
confusion_matrix <- function(score,status, cutoff = median(score), plot.it = TRUE,title = NULL) {
  
  
  
  cdscore <- ifelse(score<= cutoff, "Low", "High")
  confusion <- table(newscore,status)
  if (plot.it) fourfoldplot(confusion, color = c("#CC6666", "#99CC99"),
                            conf.level = 0, margin = 1, main = title)
  confusion
  
}

#time-dependent ROC
#get time dependent AUC dataframe of linear predictor of  markers 
getTimeDentAUCdf_multi<-function(time,status,marker,marekrstr=" ",timeline=360){
  #1 remove NA 
  df=data.frame(time,status,marker)
  df<-df[complete.cases(df),]
  a=c()
  for(i in 3:ncol(df)){
    a=c(a,paste("df[,",i,"]"))
  }
  fm<- as.formula(paste("Surv(df[,1], df[,2]) ~ ", paste(a, collapse= "+")))
  fit=coxph(fm)
  
  res=risksetROC(Stime=df[,1], status=df[,2],
                 marker=fit$linear.predictors, predict.time=timeline, method="Cox",plot = FALSE) 
  outdf<-data.frame(TP=res$TP,FP=res$FP,Marker=rep(paste(marekrstr," ",round(res$AUC,2)),length(res$FP)))
  
  return(outdf)
}
#get time dependent AUC dataframe of a indivicual  marker 
getTimeDentAUCdf<-function(time,status,marker,marekrstr=" ",timeline=360){
  
  #filter roc values
  ROCfilter<-function(x){
    #sort data by column 1
    x <-x[order(x[,1]),]
    
    #set temp maxB
    maxB <- x[1,2]
    #remove decreased B
    for (i in 1:nrow(x)){
      if(x[i,2]>=maxB){
        maxB<-x[i,2]
        #print(maxB)
      }else{
        x <- x[-i,]
        i <- i-1
        #      print(i)
      }
    }
    return(x)
  }
  
  #1 remove NA 
  df=data.frame(time,status,marker)
  df<-df[complete.cases(df),]
  #timeROC
  library(timeROC)
  df$pred <- predict(cph(Surv(df$time,df$status)~df$marker,  data=df,x=TRUE,y=TRUE,surv=TRUE))
  roc1<-timeROC(T=df$time,delta=df$status,
                marker=df$pred,cause=1,
                weighting="marginal", other_markers=NULL,
                times=timeline,iid=TRUE)
  auc<-round(roc1$AUC[2],digits = 2)  ###  AUC
  auc.ci<-round(confint(roc1)$CI_AUC/100,digits = 2) ### 95% CI of AUC
  tempdf<-data.frame(TP=roc1$TP[,2],FP=roc1$FP[,2])
  #to get concurve ROCs
  tempdf<-ROCfilter(tempdf)
  outdf<-data.frame(TP=tempdf$TP,FP=tempdf$FP,Marker=rep(paste(marekrstr," ",auc,paste(auc.ci,collapse = "-")),nrow(tempdf)))
  
  #risketROC
#   fit=coxph(Surv(df[,1], df[,2]) ~ df[,3])
#   res=risksetROC(Stime=df[,1], status=df[,2],
#                  marker=fit$linear.predictors, predict.time=timeline, method="Cox",plot = FALSE) 
#   outdf<-data.frame(TP=res$TP,FP=res$FP,Marker=rep(paste(marekrstr," AUC:",round(res$AUC,4)),length(res$FP)))
  
  list(roc1,outdf)
}

#plot time dependent ROC curves of multiple lines
getTimeDentAUCdf_matrix<-function(time,status,marker.matrix,time.line=60,data.type="DFS"){
  marker.names<-colnames(marker.matrix)
  tempdf<-data.frame(TP=c(),FP=c(),Marker=c())
  roc.list<-list()
  for (i in 1:length(marker.names) ) {
    tempdf2<-getTimeDentAUCdf(time,status,marker.matrix[,i],marekrstr=marker.names[i],timeline=time.line)
    tempdf=rbind(tempdf,tempdf2[[2]])
    roc.list[[length(roc.list)+1]]<-tempdf2[[1]]
  }
  
  #test and print pvalue between ROC curves when multivariable involved 
  if(length(marker.names)>1){
    roc1=roc.list[[1]]
    for (i in 2:length(roc.list)) {
      roc2<-roc.list[[i]]
      p <- 1-pnorm(abs(roc1$AUC[2]-roc2$AUC[2])/(((roc1$inference$vect_sd_1[2]^2+roc2$inference$vect_sd_1[2]^2)/2)^.5/2))
      print(paste(marker.names[i]," vs ",marker.names[1],": P= ", p))
    }
  }
  
  
  
  plot_time_depent_surival_AUC(tempdf,cuttime =time.line, datatype=data.type)
  
}

#plot time dependent AUC 
plot_time_depent_surival_AUC<-function(df,cuttime,datatype){
  p=ggplot(data=df,aes(x=FP,y=TP,color=Marker))+geom_line(size=1)+
    labs(title=paste("ROC for ",cuttime," months survival predicting in ",datatype,"dataset",sep=" "),x="False positive rate", y = "True negative rate")+
    scale_color_aaas()+theme_base()+theme(legend.position="top")+geom_abline(intercept = 0,linetype="dotted")+
    scale_y_continuous(limits=c(0,1),expand = c(0,0),breaks=c(0.2,0.4,0.6,0.8,1))+ 
    scale_x_continuous(limits=c(0,1),expand = c(0,0),breaks=c(0,0.2,0.4,0.6,0.8,1))+
    theme(
      legend.position = c(0.7, 0.25), # c(0,0) bottom left, c(1,1) top-right
      legend.background = element_rect(colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 12, face = 'bold')
    )
  return(p)
}


#plot ROC with multiple predictors
plotmulti_ROC_with_test<- function(predictor.matrix,status,predictor.name="predictor"){
  #loading packages 
  suppressMessages(library("ROCR"))# get colors 
  suppressMessages(library("ggsci"))# get colors
  #build data list
  pred.number<-ncol(predictor.matrix)
  col.list<-pal_npg("nrc", alpha = 1)(pred.number)
  perf.list<-list()
  legend.list<-c()
  name.list<-names(predictor.matrix)
  rocobj.list<-list()
  for(i in 1:pred.number){
    temp.vec<-predictor.matrix[,i]
    if(class(temp.vec)=="factor"){
      temp.vec<-as.numeric(temp.vec)
    }
    pred<-prediction(temp.vec,status)
    perf <- performance(pred, "tpr", "fpr" )
    auc <- round(performance(pred,"auc")@y.values[[1]][1],digits = 2)
    #make sure the relative value 
    if(auc<0.5){
      auc=1-auc
    }
    perf.list[[length(perf.list)+1]]<-perf
    #get ci 
    rocobj <- roc(status, temp.vec)
    rocobj.list[[length(rocobj.list)+1]]<-rocobj
    legend.list=c(legend.list,
                  paste(name.list[i],
                        "  ",auc,"  ",
                        paste(round(as.numeric(ci.auc(rocobj)),digits = 2)[c(1,3)],collapse="-")))
  }
  #ROC of each performance 
  ROCR::plot( perf.list[[1]], box.lty=1, box.lwd=2,col=col.list[1],
              box.col="black", lwd=2, colorkey.relwidth=0.5, xaxis.cex.axis=1,
              xaxis.col='black', xaxis.col.axis="black", yaxis.col='black', yaxis.cex.axis=1,
              yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1, xaxis.lwd=1, yaxis.lwd=1,
              ylim=c(0, 1), xlim=c(0,1),
              yaxis.col.axis="black", cex.lab=1.5, cex.main=1,title=predictor.name)
  for(i in 2:pred.number)
  {
    ROCR::plot(perf.list[[i]], add = TRUE,colorize = F,col=col.list[i],lwd=2)
  }
  print(length(legend.list))
  legend(x=0.3,y=0.3, legend = legend.list , col=col.list, lty=1, cex=1)
  abline(0,1,col="black",lty=3) 
  
  #get test pvalue 
  for(i in 2:pred.number)
  {
    p<-roc.test(rocobj.list[[1]], rocobj.list[[i]], paired=T, method="delong")
  print(paste(name.list[i]," vs ",name.list[i]," : ",round(p$p.value,digits = 4)))
  }
}
#get ROC plot of individual predictor
plotSingle_ROC_with_test<- function(predictor,status,predictor.name="predictor"){
  suppressMessages(library("ROCR"))
  pred <- prediction(predictor,status)
  pref <- performance(pred, "tpr", "fpr" )
  auc <- round(performance(pred,"auc")@y.values[[1]][1],digits = 4)
  ROCR::plot( pref, box.lty=1, box.lwd=2,col="#57719D",
              box.col="black", lwd=2, colorkey.relwidth=0.5, xaxis.cex.axis=1,
              xaxis.col='black', xaxis.col.axis="black", yaxis.col='black', yaxis.cex.axis=1,
              yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1, xaxis.lwd=1, yaxis.lwd=1,
              ylim=c(0, 1), xlim=c(0,1),
              yaxis.col.axis="black", cex.lab=1.5, cex.main=1,title=predictor.name)
  legend(x=0.5,y=0.2, legend =paste(predictor.name," AUC =",auc) , col='#57719D', lty=1, cex=1)
  abline(0,1,col="red") 
}

#get print version of HR 
get_print_HR_cox<-function(surv,variable,rev=F){
  if(rev){
    variable<-factor(variable,level=rev(levels(as.factor(variable))))
  }
  df<-prettify(summary(coxph(surv~variable)))
  
  a=paste(round(df$`Hazard Ratio`,digits = 2)," (",round(df$`CI (lower)`,digits = 2),"-",round(df$`CI (upper)`,digits = 2),")",sep="")
  return(a)
  
}


#customized survival plot 

cust_surp<-function(fit,df){
  ggsurvplot(fit, data = df,
             surv.median.line = "hv", # Add medians survival
             # Change legends: title & labels
             legend.title = "Risk",
             legend.labs = c("Low", "High"),
             # Add p-value and confidence intervals
             pval = TRUE,
             conf.int = TRUE,
             # Add risk table
             risk.table = TRUE,
             tables.height = 0.2,
             tables.theme = theme_cleantable(),
             # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
             # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
             palette = c("#E7B800", "#2E9FDF"),
             ggtheme = theme_bw() # Change ggplot2 theme
  )
}


Initial_surivival <- function(dataWithCli,type="DFS",plot.title=""){
  
  .formula <<- as.formula("Surv(dataWithCli$DFS, dataWithCli$DFSstatus==1)~ dataWithCli$riskscoreStatus")
  if(type=="OS"){
    .formula <- as.formula("Surv(dataWithCli$OS, dataWithCli$OSstatus==1)~ dataWithCli$riskscoreStatus")
  }
  surv.fit<- survfit(.formula)
  survtotal.diff<- survdiff(.formula)
  survtotal.pvalue<-1 - pchisq(survtotal.diff$chisq, length(survtotal.diff$n) - 1)
  legend.position=c(0.8,0.3)
  surpt.pt.dfs=ggsurvplot(
    surv.fit,
    pval = T,
    break.time.by = 12,
    risk.table = T,
#     risk.table.fontsize = 6,
#     tables.height = 0.15,
#     surv.plot.height=3,
    tables.theme = theme_survminer( font.legend  = c(3, "plain", "darkgreen")),
    # risk.table.pos="in",
    # Useful when you have multiple groups
    censor = T,
    legend.title = "Riskscore",
    legend.labs = c("High", "Low"),
    palette = c("#a50026","blue"),
    legend=legend.position,
    # 
    title=plot.title,
    ggtheme = theme_classic()
  )
  return(surpt.pt.dfs)
}

quote_var<-function(value){
  for (i in 1:length(value)) {
    value[i] <- paste("`",value[i],"`",sep = "")
  }
  return(value)
}

cust_surp2<-function(fit,df=NULL,colorlist=NULL,ti=""){
  if(is.null(colorlist)){
    ggsurvplot(fit, data = df,
               surv.median.line = "hv", # Add medians survival
               
               # Change legends: title & labels
               legend.title = "Risk",
               # Add p-value and confidence intervals
               pval = TRUE,
               title=ti,
               legend.labs = c("Low", "High"),
               # Add p-value and confidence intervals
               conf.int = F,
               # Add risk table
               risk.table = TRUE,
               tables.height = 0.2,
               tables.theme = theme_cleantable(),
               palette = c("#E7B800", "#2E9FDF"),
               # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
               # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
               ggtheme = theme_bw() # Change ggplot2 theme
    )
  }else{
    ggsurvplot(fit, data = df,
               surv.median.line = "hv", # Add medians survival
               
               # Change legends: title & labels
               legend.title = "Risk",
               # Add p-value and confidence intervals
               legend.labs = c("Low", "High"),
               # Add p-value and confidence intervals
               conf.int = F,
               title=ti,
               # Add risk table
               risk.table = TRUE,
               tables.height = 0.2,
               tables.theme = theme_cleantable(),
               palette = c("#E7B800", "#2E9FDF"),               # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
               # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
               ggtheme = theme_bw() # Change ggplot2 theme
    )
  }
  
}

#plot time-dependent ROC 



