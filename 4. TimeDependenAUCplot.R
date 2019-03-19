library("survivalROC")
library("survival")
library("RColorBrewer")
library(survminer)

get_liner_predictor<-function(combineMaker,time.dependent.data){
  formula<-as.formula(paste("Surv(time.dependent.data$DFS,time.dependent.data$DFSstatus)~",paste(combineMaker,collapse = "+"),sep=""))
  coxfit<-coxph(formula,data=time.dependent.data)
  time.dependent.data.new<-cbind(time.dependent.data,coxfit$linear.predictors)
  auc.temp<-c()
  for(t in 6:120){
    auc.temp=c(auc.temp,survivalROC.C(Stime=time.dependent.data$OS,status =time.dependent.data$OSstatus,marker =coxfit$linear.predictors,predict.time = t,span = 0.05)$AUC)
  }
  return(auc.temp)
}


names(trainingData)
time.dependent.data<-ValidationData[,c(lassoMakers_DFS,"riskscore","DFS","DFSstatus","OS","OSstatus")]
auc.1=NULL
auc.2=NULL
auc.3=NULL
auc.4=NULL
auc.cirscore=NULL
sur.dept<-survival::Surv(time.dependent.data$DFS,time.dependent.data$DFSstatus==1)
for(t in 6:120){
  auc.1=c(auc.1,survivalROC.C(Stime=time.dependent.data$DFS,status =time.dependent.data$DFSstatus,marker =time.dependent.data[,1],predict.time = t,span = 0.05)$AUC)
  auc.2=c(auc.2,survivalROC.C(Stime=time.dependent.data$DFS,status =time.dependent.data$DFSstatus,marker =time.dependent.data[,2],predict.time = t,span = 0.05)$AUC)
  auc.3=c(auc.3,survivalROC.C(Stime=time.dependent.data$DFS,status =time.dependent.data$DFSstatus,marker =time.dependent.data[,3],predict.time = t,span = 0.05)$AUC)
  auc.4=c(auc.4,survivalROC.C(Stime=time.dependent.data$DFS,status =time.dependent.data$DFSstatus,marker =time.dependent.data[,4],predict.time = t,span = 0.05)$AUC)
  auc.cirscore=c(auc.cirscore,survivalROC.C(Stime=time.dependent.data$DFS,status =time.dependent.data$DFSstatus,marker =time.dependent.data[,5],predict.time = t,span = 0.05)$AUC)
  
}
#plot 
plot(6:120,auc.1,type="l",xlab="Time",ylab="AUC",col="red",ylim = c(0.5,1),lty=2,lwd=2)
lines(6:120,1-auc.2,col="blue",lty=2,lwd=2)
lines(6:120,auc.3,col="green",lty=2,lwd=2)
lines(6:120,auc.4,col="purple",lty=2,lwd=2)
lines(6:120,auc.cirscore,col="black",lty=1,lwd=2)
legend(80,1,legend = c("cirScore",lassoMakers_DFS),col=c("black","red","blue","green","purple"),lty=c(1,2,2,2,2),bt="n",lwd=2)

# test the significance between cirScore and individual scores
wilcox.test(auc.cirscore,auc.1)
wilcox.test(auc.cirscore,auc.2)
wilcox.test(auc.cirscore,auc.3)
wilcox.test(auc.cirscore,auc.4)


# test for other combination


comb.2marker<-combn(lassoMakers_DFS,2)
comb.3marker<-combn(lassoMakers_DFS,3)
color.2marker<-brewer.pal(ncol(comb.2marker), "Set1")
color.3marker<-brewer.pal(ncol(comb.3marker), "Set2")

# all combine data 
auc.cirscore=get_liner_predictor(lassoMakers_DFS,time.dependent.data)


# plot two markers
auc.list2<-as.data.frame(apply(comb.2marker, 2, get_liner_predictor,time.dependent.data))
auc.list2.names<-apply(comb.2marker, 2, paste,collapse="/")
colnames(auc.list2)<-auc.list2.names
plot(6:120,auc.cirscore,type="l",xlab="Time",ylab="AUC",col="black",ylim = c(0.5,1),lty=1,lwd=2)
for (i in 1:ncol(auc.list2)) {
  lines(6:120,auc.list2[,i],col=color.2marker[i],lty=2,lwd=2)
}
legend(20,1,legend = c("cirScore",auc.list2.names),col=c("black",color.2marker),lty=c(1,rep(2,ncol(auc.list2))),bt="n",lwd=2)

# plot three markers
auc.list3<-as.data.frame(apply(comb.3marker, 2, get_liner_predictor,time.dependent.data))
auc.list3.names<-apply(comb.3marker, 2, paste,collapse="/")
colnames(auc.list3)<-auc.list3.names
plot(6:120,auc.cirscore,type="l",xlab="Time",ylab="AUC",col="black",ylim = c(0.5,1),lty=1,lwd=1)
for (i in 1:ncol(auc.list3)) {
  lines(6:120,auc.list3[,i],col=color.3marker[i],lty=1,lwd=1)
}
legend(10,1,legend = c("cirScore",auc.list3.names),col=c("black",color.3marker),lty=c(1,rep(1,ncol(auc.list2))),bt="n",lwd=2)

# ggforest 



formula<-as.formula(paste("Surv(trainingData$DFS,trainingData$DFSstatus)~",paste(lassoMakers_DFS,collapse = "+"),sep=""))
coxfit<-coxph(formula,data=trainingData)
ggforest(coxfit,data = trainingData)


