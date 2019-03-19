# combine training data and internal validation data 
stratify.variable<-c("riskscoreStatus","stage_risk","OS","OSstatus","DFS","DFSstatus")
#stage 
stratifi.combinedata<-rbind(trainingData[,stratify.variable],ValidationData[,stratify.variable],extervaliData[,stratify.variable])
stratifi.combinedata$stage_risk<-as.factor(stratifi.combinedata$stage_risk)
stratifi.combinedata$riskscoreStatus=factor(stratifi.combinedata$riskscoreStatus,levels=rev(levels(as.factor(stratifi.combinedata$riskscoreStatus))))
stratifi.combinedata$combineF<-with(stratifi.combinedata,interaction(stage_risk,riskscoreStatus))
vec.res<-levels(stratifi.combinedata$stage_risk)
#DFS
dfs.list<-list()
for(i in 1:length(vec.res)){
  sub.data<-stratifi.combinedata[stratifi.combinedata$stage_risk==vec.res[i],]
  subgroup.fit.dfs<-survfit(Surv(sub.data$DFS,sub.data$DFSstatus==1)~ riskscoreStatus,sub.data)
  pdf(paste("20180626/strat.",vec.res[i],".dfs.pdf",sep=""),width=8,height=8) 
  tempp<-cust_surp2(subgroup.fit.dfs,sub.data)
  print(tempp,newpage=FALSE)
  dev.off()
  print(get_print_HR_cox(Surv(sub.data$DFS,sub.data$DFSstatus==1),sub.data$riskscoreStatus))
  print(vec.res[i])
  
  # dfs.list[[length(dfs.list)+1]] <- subp
}
#OS
os.list<-list()
for(i in 1:length(vec.res)){
  sub.data<-stratifi.combinedata[stratifi.combinedata$stage_risk==vec.res[i],]
  subgroup.fit.dfs<-survfit(Surv(sub.data$OS,sub.data$OSstatus==1)~ riskscoreStatus,sub.data)
  pdf(paste("20180626/strat.",vec.res[i],".os.pdf",sep=""),width=8,height=8) 
  tempp<-cust_surp2(subgroup.fit.dfs,sub.data)
  print(tempp,newpage=FALSE)
  dev.off()
  print(get_print_HR_cox(Surv(sub.data$OS,sub.data$OSstatus==1),sub.data$riskscoreStatus))
  print(vec.res[i])
  
  # os.list[[length(dfs.list)+1]] <- subp
}
# subgroup.fit.dfs<-coxph(Surv(stratifi.combinedata$DFS,stratifi.combinedata$DFSstatus==1) ~riskscoreStatus+strata(stage_risk),data = stratifi.combinedata )


#refactor levels 
subgroup.fit.dfs<-survfit(Surv(stratifi.combinedata$DFS,stratifi.combinedata$DFSstatus==1)~ combineF,stratifi.combinedata)
p.sub.1=cust_surp(subgroup.fit.dfs,stratifi.combinedata)
#show hr
subgroup.fit.os<-survfit(stratifi.combinedata.os~ riskscoreStatus,stratifi.combinedata)
p.sub.2=cust_surp(subgroup.fit.os,stratifi.combinedata)
print("Stage II")
get_print_HR_cox(stratifi.combinedata.dfs,stratifi.combinedata$riskscoreStatus)
get_print_HR_cox(stratifi.combinedata.os,stratifi.combinedata$riskscoreStatus)

#stage III
stratifi.combinedata<-rbind(trainingData,ValidationData)
stratifi.combinedata<-stratifi.combinedata[stratifi.combinedata$Stage==3,]
stratifi.combinedata$riskscoreStatus<-ifelse(stratifi.combinedata$riskscore>=cutoff,"high","low")
stratifi.combinedata$riskscoreStatus=factor(stratifi.combinedata$riskscoreStatus,levels=rev(levels(as.factor(stratifi.combinedata$riskscoreStatus))))
stratifi.combinedata.dfs<-Surv(stratifi.combinedata$DFS,stratifi.combinedata$DFSstatus)
stratifi.combinedata.os<-Surv(stratifi.combinedata$OS,stratifi.combinedata$OSstatus)
subgroup.fit.dfs<-survfit(stratifi.combinedata.dfs~ riskscoreStatus,stratifi.combinedata)
p.sub.3=cust_surp(subgroup.fit.dfs,stratifi.combinedata)

subgroup.fit.os<-survfit(stratifi.combinedata.os~ riskscoreStatus,stratifi.combinedata)
p.sub.4=cust_surp(subgroup.fit.os,stratifi.combinedata)
print("Stage III")
get_print_HR_cox(stratifi.combinedata.dfs,stratifi.combinedata$riskscoreStatus)
get_print_HR_cox(stratifi.combinedata.os,stratifi.combinedata$riskscoreStatus)


#
savePlot(p.sub.1)
savePlot(p.sub.2)
savePlot(p.sub.3)
savePlot(p.sub.4)
plot_grid(plot_grid(p.sub.1$plot,p.sub.1$table,ncol=1),plot_grid(p.sub.2$plot,p.sub.2$table,ncol=1)
          ,plot_grid(p.sub.3$plot,p.sub.3$table,ncol=1),plot_grid(p.sub.4$plot,p.sub.4$table,ncol=1),
          ncol = 2, labels = c("A", "B","C","D"))
#get HR


