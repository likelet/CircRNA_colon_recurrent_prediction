#preloading functions and packages ------
source("./analysis_functions.R")

# Function code for scale data and permutation 
scale_impute_data<-function(expdata){
  tempdata=expdata
  tempdata=log2(tempdata*10000)
  #identifi outlie 
  qnt <- quantile(tempdata, probs=c(.25, .75), na.rm = T)
  H <- 1.5 * IQR(as.matrix(tempdata), na.rm = T)
  tempdata=as.matrix(tempdata)
  tempdata[tempdata < (qnt[1] - H) | tempdata >(qnt[2] + H)] <- NA
  #do permutation 
  tempdata.impute=impute.knn(as.matrix(tempdata) ,k = 10, rowmax = 0.5, colmax = 0.8)
  normalize.Data=as.data.frame(tempdata.impute$data)
  return(normalize.Data)
}





trainingData=read.csv("data/trainingData.csv",row.names=1,header=T)

#DFS for unicox regression analysis 
resultTableCox=batchUnivarCOXfun(Surv(round(trainingData$OS), trainingData$OSstatus==1),trainingData[,2:22])
marker_candidate=row.names(resultTableCox)[resultTableCox[,7]<0.05]

#heatmap functions
pheatmap(trainingData[,marker_candidate],scale="none",fontsize_row = 5) 

#DFS
selecVlist1=c()
for(i in 1:500){
  sampleindex2=sample(1:nrow(trainingData),1*nrow(trainingData),rep=T)
  effectdata=datawithCli[sampleindex2,]
  
  glmmod<-glmnet(as.matrix(effectdata[,marker_candidate]),Surv(round(effectdata$DFS), effectdata$DFSstatus==1), family = "cox")
  cv.glmmod<-cv.glmnet(as.matrix(effectdata[,marker_candidate]),Surv(round(effectdata$DFS), effectdata$DFSstatus==1), family = "cox")
  best_lambda <- cv.glmmod$lambda.1se
  result<-coef(glmmod, s = best_lambda)
  selecVlist1=c(selecVlist1,result@Dimnames[[1]][which(result != 0)])
  print(i)
}
tablecount1=table(selecVlist1)
par(mar=c(9,4,2,2))
# 
# #plot frequent table
tablecountdf=data.frame(rev(sort(tablecount1))/500)
names(tablecountdf)[1]="GeneName"
plottable=ggplot(tablecountdf,aes(x = GeneName, y = Freq))+
  geom_col() +
  ylab ("Frequency of markers selected")+geom_hline(yintercept=0.2,col="red")+ 
  scale_y_continuous(limit=c(0,1),expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.1,size=10))
#get lasso names 
lassoMakers_DFS <- names(tablecount1[tablecount1>=75])

fmla <- as.formula(paste("Surv(trainingData$DFS, trainingData$DFSstatus==1) ~ ", paste(lassoMakers_DFS, collapse= "+")))
dfs.cox.fit2 <- coxph(fmla, data=trainingData)
trainingData$riskscore=predict(dfs.cox.fit2,trainingData)

# similarly, we calculated the riskscore(cirScore) for internal validation data and external validation dataset 







