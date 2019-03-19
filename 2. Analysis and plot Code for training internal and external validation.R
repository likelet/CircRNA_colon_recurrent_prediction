



#preloading function------
source("./analysis_functions.R")

trainingData <- fread("data/trainingWithCirScore.csv",header = T,check.names = F)
ValidationData <- fread("data/intvalidWithCirScore.csv",header = T,check.names = F)
extervaliData <- fread("data/exvalid_180WithCirScore.csv",header = T,check.names = F)
#catogrized data set ----
    trainingData$ageStatus <- ifelse(trainingData$Age>=65,1,0)
    ValidationData$ageStatus <- ifelse(ValidationData$Age>=65,1,0)
    extervaliData$ageStatus <- ifelse(extervaliData$Age>=65,1,0)
    #Positive Node
    trainingData$pN_status<-I(trainingData$pN>=1)+I(trainingData$pN>=4)
    ValidationData$pN_status<-I(ValidationData$pN>=1)+I(ValidationData$pN>=4)
    extervaliData$pN_status<-I(extervaliData$pN>=1)+I(extervaliData$pN>=4)
    #total Node
    trainingData$Ntotal_status<-ifelse(trainingData$Ntotal>=12,1,0)
    ValidationData$Ntotal_status<-ifelse(ValidationData$Ntotal>=12,1,0)
    extervaliData$Ntotal_status<-ifelse(extervaliData$Ntotal>=12,1,0)
    #T_stage 
    trainingData$Tstage_status<-ifelse(trainingData$T_stage>=4,1,0)
    ValidationData$Tstage_status<-ifelse(ValidationData$T_stage>=4,1,0)
    extervaliData$Tstage_status<-ifelse(extervaliData$T_stage=="T4",1,0)
    #differential 
    trainingData$grade_status<-ifelse(trainingData$grade>=3,1,0)
    ValidationData$grade_status<-ifelse(ValidationData$grade>=3,1,0)
    extervaliData$grade_status<-ifelse(extervaliData$grade>=3,1,0)

    trainingData<-data.frame(trainingData,check.names = F)
    ValidationData<-data.frame(ValidationData,check.names = F)
    extervaliData<-data.frame(extervaliData,check.names = F)
    #set surv object 
    survtotal.train.dfs <- Surv(trainingData$DFS, trainingData$DFSstatus==1)
    survtotal.train.os <- Surv(trainingData$OS, trainingData$OSstatus==1)
    survtotal.invali.dfs <- Surv(ValidationData$DFS, ValidationData$DFSstatus==1)
    survtotal.invali.os <- Surv(ValidationData$OS, ValidationData$OSstatus==1)
    survtotal.exvali.dfs <- Surv(extervaliData$DFS, extervaliData$DFSstatus==1)
    survtotal.exvali.os <- Surv(extervaliData$OS, extervaliData$OSstatus==1)
#plot survival curve of cirScore ---- 
        library("survminer")
        
        #KM curves 
        surpt.train.dfs=Initial_surivival(trainingData,type ="DFS",plot.title = "Training DFS")
        surpt.train.os=Initial_surivival(trainingData,type ="OS",plot.title = "Training OS")
        pdf("result/vali.dfs.pdf",width=8,height=8)
        Initial_surivival(ValidationData,type ="DFS",plot.title = "Internal validation DFS")
        dev.off()
        pdf("result/vali.os.pdf",width=8,height=8)
        Initial_surivival(ValidationData,type ="OS",plot.title = "Internal validation  OS")
        dev.off()
        pdf("result/exvali.dfs.pdf",width=8,height=8)
        Initial_surivival(extervaliData,type ="DFS",plot.title = "External validation DFS")
        dev.off()
        pdf("result/exvali.os.pdf",width=8,height=8)
        Initial_surivival(extervaliData,type ="OS",plot.title = "External validation OS")
        dev.off()
        

        
        
        #get HR ratio
        get_print_HR_cox(survtotal.train.dfs,trainingData$riskscoreStatus,rev=T)
        get_print_HR_cox(survtotal.train.os,trainingData$riskscoreStatus,rev=T)
        get_print_HR_cox(survtotal.invali.dfs,ValidationData$riskscoreStatus,rev=T)
        get_print_HR_cox(survtotal.invali.os,ValidationData$riskscoreStatus,rev=T)
        get_print_HR_cox(survtotal.exvali.dfs,extervaliData$riskscoreStatus,rev=T)
        get_print_HR_cox(survtotal.exvali.os,extervaliData$riskscoreStatus,rev=T)
        



# multi variable analysis----
        interestVariable<-c("Sex","ageStatus","LeftOrRight","NI","VI","pN_status","Ntotal_status","Tstage_status","grade_status","riskscoreStatus")
        
        trainingData$ageStatus=as.factor(trainingData$ageStatus)
        trainingData$NI=as.factor(trainingData$NI)
        trainingData$VI=as.factor(trainingData$VI)
        trainingData$Stage=as.factor(trainingData$Stage)
        trainingData$riskscoreStatus=as.factor(trainingData$riskscoreStatus)
        
        #unicox 
        unicoxtable=batchUnivarCOXfun(Surv(trainingData$DFS, trainingData$DFSstatus==1),trainingData[,which(names(trainingData) %in% interestVariable)])
        unicoxtable=round(unicoxtable,digits = 3)
        marker_unicox_cli=row.names(unicoxtable)[unicoxtable[,7]<0.05]
        # write.csv(unicoxtable,paste("data/",Sys.time(),"training_DFS_unicox.csv",sep = ""))
        #multivarialbe cox regression
        #trainingDFS
        trainingData$riskscoreStatus=factor(trainingData$riskscoreStatus,levels = rev(levels(trainingData$riskscoreStatus)))
        multifmla.train.DFS <- as.formula(paste("Surv(trainingData$DFS, trainingData$DFSstatus==1) ~ ", paste(interestVariable, collapse= "+")))
        cox.fit2_multi.train.DFS <- coxph(multifmla.train.DFS , data=trainingData)
        cox.fit2_multi.train.DFS.table=prettify(summary(cox.fit2_multi.train.DFS),digits = 3)
        cox.fit2_multi.train.DFS.table.fw=zhaoqi.write.csv.survival(cox.fit2_multi.train.DFS.table)
        cox.fit2_multi.train.DFS.table.fw
        write.csv(cox.fit2_multi.train.DFS.table.fw,paste("result/","training_DFS_multicox.csv",sep = ""))
        #traningOS
        multifmla.train.OS <- as.formula(paste("Surv(trainingData$OS, trainingData$OSstatus==1) ~ ", paste(interestVariable, collapse= "+")))
        cox.fit2_multi.train.OS <- coxph(multifmla.train.OS, data=trainingData)
        cox.fit2_multi.train.OS.table=prettify(summary(cox.fit2_multi.train.OS),digits = 3)
        cox.fit2_multi.train.OS.table.fw=zhaoqi.write.csv.survival(cox.fit2_multi.train.OS.table)
        cox.fit2_multi.train.OS.table.fw
        write.csv(cox.fit2_multi.train.OS.table.fw,paste("result/","training_OS_multicox.csv",sep = ""))
        #valiDFS
        ValidationData$ageStatus=as.factor(ValidationData$ageStatus)
        ValidationData$NI=as.factor(ValidationData$NI)
        ValidationData$VI=as.factor(ValidationData$VI)
        ValidationData$Stage=as.factor(ValidationData$Stage)
        ValidationData$riskscoreStatus=as.factor(ValidationData$riskscoreStatus)
        ValidationData$riskscoreStatus=factor(ValidationData$riskscoreStatus,levels = rev(levels(ValidationData$riskscoreStatus)))
        multifmla.vali.DFS <- as.formula(paste("Surv(ValidationData$DFS, ValidationData$DFSstatus==1) ~ ", paste(interestVariable, collapse= "+")))
        
        cox.fit2_multi.vali.DFS <- coxph(multifmla.vali.DFS, data=ValidationData)
        cox.fit2_multi.vali.DFS.table=prettify(summary(cox.fit2_multi.vali.DFS),digits = 3)
        cox.fit2_multi.vali.DFS.table.fw=zhaoqi.write.csv.survival(cox.fit2_multi.vali.DFS.table)
        cox.fit2_multi.vali.DFS.table.fw
        write.csv(cox.fit2_multi.vali.DFS.table.fw,paste("result/","vali_DFS_multicox.csv",sep = ""))
        #valiOS
        multifmla.vali.OS <- as.formula(paste("Surv(ValidationData$OS, ValidationData$OSstatus==1) ~ ", paste(quote_var(interestVariable), collapse= "+")))
        cox.fit2_multi.vali.OS <- coxph(multifmla.vali.OS, data=ValidationData)
        cox.fit2_multi.train.vali.table=prettify(summary(cox.fit2_multi.vali.OS),digits = 3)
        cox.fit2_multi.train.vali.table.forw=zhaoqi.write.csv.survival(cox.fit2_multi.train.vali.table)
        cox.fit2_multi.train.vali.table.forw
        write.csv(cox.fit2_multi.train.vali.table.forw,paste("result/","vali_OS_multicox.csv",sep = ""))
        
        #extvaliDFS
        extervaliData$ageStatus=as.factor(extervaliData$ageStatus)
        extervaliData$NI=as.factor(extervaliData$NI)
        extervaliData$VI=as.factor(extervaliData$VI)
        extervaliData$Stage=as.factor(extervaliData$Stage)
        extervaliData$riskscoreStatus=as.factor(extervaliData$riskscoreStatus)
        extervaliData$riskscoreStatus=factor(extervaliData$riskscoreStatus,levels = rev(levels(extervaliData$riskscoreStatus)))
        multifmla.exvali.DFS <- as.formula(paste("Surv(extervaliData$DFS, extervaliData$DFSstatus==1) ~ ", paste(interestVariable, collapse= "+")))
        cox.fit2_multi.exvali.DFS <- coxph(multifmla.exvali.DFS, data=extervaliData)
        cox.fit2_multi.exvali.DFS.table=prettify(summary(cox.fit2_multi.exvali.DFS),digits = 3)
        cox.fit2_multi.exvali.DFS.table.fw=zhaoqi.write.csv.survival(cox.fit2_multi.exvali.DFS.table)
        cox.fit2_multi.exvali.DFS.table.fw
        write.csv(cox.fit2_multi.exvali.DFS.table.fw,paste("result/","exvali_DFS_multicox.csv",sep = ""))
        #exvaliOS
        multifmla.exvali.OS <- as.formula(paste("Surv(extervaliData$OS, extervaliData$OSstatus==1) ~ ", paste(interestVariable, collapse= "+")))
        cox.fit2_multi.exvali.OS <- coxph(multifmla.exvali.OS, data=extervaliData)
        cox.fit2_multi.train.exvali.table=prettify(summary(cox.fit2_multi.exvali.OS),digits = 3)
        cox.fit2_multi.train.exvali.table.forw=zhaoqi.write.csv.survival(cox.fit2_multi.train.exvali.table)
        cox.fit2_multi.train.exvali.table.forw
        write.csv(cox.fit2_multi.train.exvali.table.forw,paste("result/","exvali_OS_multicox.csv",sep = ""))
#----------
########ROC analysis --------
        
        #roc analysis for individual data 
        #traing data 
        interestVariable.ROC<-c("riskscore","","","NI","VI")
        trainingData<-data.frame(trainingData)
        traindataforROC<-trainingData[,interestVariable.ROC]
        traindataforROC$Sex=ifelse(traindataforROC$Sex=="Male",1,0)
        pdf("result/ROC.train.dfs.pdf",width=8,height=8) 
        plotmulti_ROC_with_test(traindataforROC,trainingData$DFSstatus)
        dev.off()
        pdf("result/ROC.train.os.pdf",width=8,height=8) 
        plotmulti_ROC_with_test(traindataforROC,trainingData$OSstatus)
        dev.off()
        #internal validation 
        ValidationData<-as.data.frame(ValidationData)
        vali.forROC<-ValidationData[,interestVariable.ROC]
        vali.forROC$Sex=ifelse(vali.forROC$Sex=="Male",1,0)
        pdf("result/ROC.vali.dfs.pdf",width=8,height=8) 
        plotmulti_ROC_with_test(vali.forROC,ValidationData$DFSstatus,predictor.name = "Vali.DFS")
        dev.off()
        pdf("result/ROC.vali.os.pdf",width=8,height=8) 
        plotmulti_ROC_with_test(vali.forROC,ValidationData$OSstatus,predictor.name = "Vali.OS")
        dev.off()
        # external validation
        extervaliData<-as.data.frame(extervaliData)
        external.vali.forROC<-extervaliData[,interestVariable.ROC]
        external.vali.forROC$Sex=ifelse(external.vali.forROC$Sex=="Male",1,0)
        pdf("result/ROC.exvali.dfs.pdf",width=8,height=8) 
        plotmulti_ROC_with_test(external.vali.forROC,extervaliData$DFSstatus,predictor.name = "exVali.DFS")
        dev.off()
        pdf("result/ROC.exvali.os.pdf",width=8,height=8) 
        plotmulti_ROC_with_test(external.vali.forROC,extervaliData$OSstatus,predictor.name = "exVali.OS")    
        dev.off()
 #--------       
#cirScore+23risk+combine ROC ----
        trainingData$stage_risk.fac<-as.factor(trainingData$stage_risk)
        ValidationData$stage_risk.fac<-as.factor(ValidationData$stage_risk)
        extervaliData$stage_risk.fac<-as.factor(extervaliData$stage_risk)
        
       
        #predict fil
        trainingData$cir.23risk.combine<-predict(cph(survtotal.train.dfs~catg(stage_risk) +catg(riskscoreStatus),trainingData,x=TRUE,y=TRUE,surv=TRUE))
        ValidationData$cir.23risk.combine<-predict(cph(survtotal.invali.dfs~catg(stage_risk) +catg(riskscoreStatus),ValidationData,x=TRUE,y=TRUE,surv=TRUE))
        extervaliData$cir.23risk.combine<-predict(cph(survtotal.exvali.dfs~catg(stage_risk) +catg(riskscoreStatus),extervaliData,x=TRUE,y=TRUE,surv=TRUE))
        interestVariable.ROC<-c("stage_risk","cir.23risk.combine")
        #get ROC data matrix 
        p1=getTimeDentAUCdf_matrix(trainingData$DFS,trainingData$DFSstatus,trainingData[,interestVariable.ROC],time.line=60,data.type = "train.DFS")
        p2=getTimeDentAUCdf_matrix(ValidationData$DFS,ValidationData$DFSstatus,ValidationData[,interestVariable.ROC],time.line=60,data.type = "vali.DFS")
        p3=getTimeDentAUCdf_matrix(extervaliData$DFS,extervaliData$DFSstatus,extervaliData[,interestVariable.ROC],time.line=36,data.type = "exvali.DFS")
        #OS prediction 
        trainingData$cir.23risk.combine.os<-predict(cph(survtotal.train.os~catg(stage_risk) +catg(riskscoreStatus),trainingData,x=TRUE,y=TRUE,surv=TRUE))
        ValidationData$cir.23risk.combine.os<-predict(cph(survtotal.invali.os~catg(stage_risk) +catg(riskscoreStatus),ValidationData,x=TRUE,y=TRUE,surv=TRUE))
        extervaliData$cir.23risk.combine.os<-predict(cph(survtotal.exvali.os~catg(stage_risk) +catg(riskscoreStatus),extervaliData,x=TRUE,y=TRUE,surv=TRUE))
        interestVariable.ROC.os<-c("stage_risk","cir.23risk.combine.os")
        #get ROC data matrix 
        p4=getTimeDentAUCdf_matrix(trainingData$OS,trainingData$OSstatus,trainingData[,interestVariable.ROC.os],time.line=60,data.type = "train.OS")
        p5=getTimeDentAUCdf_matrix(ValidationData$OS,ValidationData$OSstatus,ValidationData[,interestVariable.ROC.os],time.line=60,data.type = "vali.OS")
        p6=getTimeDentAUCdf_matrix(extervaliData$OS,extervaliData$OSstatus,extervaliData[,interestVariable.ROC.os],time.line=36,data.type = "exvali.OS")
        
        pdf("result/ROCcomparison2.pdf",width=8.27,height=11.69,paper='a4',pointsize=2) 
        plot_grid(p1,p4,p2,p5,p3,p6,ncol = 2, labels = c("A", "B","C","D","E","F"))
        dev.off()
        
        
        #combine all clinical features+individual only----
        #DFS
        interest.variable.nomo.os<-c("nomo.combine.dfs","ageStatus","NI","VI","pN_status","cirScore.status")
        trainingData$cirScore.status<-ifelse(trainingData$riskscoreStatus=="High",1,0)
        ValidationData$cirScore.status<-ifelse(ValidationData$riskscoreStatus=="High",1,0)
        extervaliData$cirScore.status<-ifelse(extervaliData$riskscoreStatus=="High",1,0)
        
        trainingData$nomo.combine.dfs<-predict(cph(survtotal.train.dfs~rev(catg(ageStatus)) +catg(NI)+catg(VI)+catg(pN_status)+catg(cirScore.status),
                                                    trainingData,x=TRUE,y=TRUE,surv=TRUE))
        ValidationData$nomo.combine.dfs<-predict(cph(survtotal.invali.dfs~catg(ageStatus) +catg(NI)+catg(VI)+catg(pN_status)+catg(cirScore.status),
                                                     ValidationData,x=TRUE,y=TRUE,surv=TRUE))
        extervaliData$nomo.combine.dfs<-predict(cph(survtotal.exvali.dfs~catg(ageStatus) +catg(NI)+catg(VI)+catg(pN_status)+catg(cirScore.status),
                                                    extervaliData,x=TRUE,y=TRUE,surv=TRUE))
        
        p1.nomo=getTimeDentAUCdf_matrix(trainingData$DFS,trainingData$DFSstatus,trainingData[,interest.variable.nomo.os],time.line=60,data.type = "train.DFS")
        p2.nomo=getTimeDentAUCdf_matrix(ValidationData$DFS,ValidationData$DFSstatus,ValidationData[,interest.variable.nomo.os],time.line=60,data.type = "vali.DFS")
        p3.nomo=getTimeDentAUCdf_matrix(extervaliData$DFS,extervaliData$DFSstatus,extervaliData[,interest.variable.nomo.os],time.line=60,data.type = "exvali.DFS")
        #OS
        interest.variable.nomo.os<-c("nomo.combine.os","ageStatus","NI","VI","pN_status","cirScore.status")
        trainingData$nomo.combine.os<-predict(cph(survtotal.train.os~catg(ageStatus) +catg(NI)+catg(VI)+catg(pN_status)+catg(cirScore.status),
                                                   trainingData,x=TRUE,y=TRUE,surv=TRUE))
        ValidationData$nomo.combine.os<-predict(cph(survtotal.invali.os~catg(ageStatus) +catg(NI)+catg(VI)+catg(pN_status)+catg(cirScore.status),
                                                    ValidationData,x=TRUE,y=TRUE,surv=TRUE))
        extervaliData$nomo.combine.os<-predict(cph(survtotal.exvali.os~catg(ageStatus) +catg(NI)+catg(VI)+catg(pN_status)+catg(cirScore.status),
                                                    extervaliData,x=TRUE,y=TRUE,surv=TRUE))
        p4.nomo=getTimeDentAUCdf_matrix(trainingData$OS,trainingData$OSstatus,trainingData[,interest.variable.nomo.os],time.line=60,data.type = "train.OS")
        p5.nomo=getTimeDentAUCdf_matrix(ValidationData$OS,ValidationData$OSstatus,ValidationData[,interest.variable.nomo.os],time.line=60,data.type = "vali.OS")
        p6.nomo=getTimeDentAUCdf_matrix(extervaliData$OS,extervaliData$OSstatus,extervaliData[,interest.variable.nomo.os],time.line=36,data.type = "exvali.OS")
        
        plot_grid(p1.nomo,p4.nomo,p2.nomo,p5.nomo,p3.nomo,p6.nomo,ncol = 2, labels = c("A", "B","C","D","E","F"))
        
        
#nomogram -------
#nomogram for DFS
        library(rms)
        data.for.stepwise<-trainingData[,c(interestVariable,"OS","OSstatus")]
        
        #backward variable selection 
        formular.stepwise <- as.formula(paste("Surv(data.for.stepwise$OS, data.for.stepwise$OSstatus==1) ~ ", paste(interestVariable, collapse= "+")))
        step(coxph(formular.stepwise,data=data.for.stepwise),direction="backward")
        step.variable<-c("ageStatus","NI","VI","pN_status","riskscoreStatus","DFS","DFSstatus")
        #package data 
        data.for.nomogram<-data.for.stepwise[,step.variable]
        dd <- datadist(data.for.nomogram)
        options(datadist="dd")
        #regression data 
        f <- cph(Surv(data.for.nomogram$DFS, data.for.nomogram$DFSstatus) ~ riskscoreStatus+ageStatus+pN_status+NI+VI , x=T, y=T, surv=T, data=data.for.nomogram, time.inc=60)
        surv <- Survival(f)
        #nomogram
        nom <- nomogram(f, fun=list(function(x) surv(36, x),function(x) surv(60, x)),  
                        lp=F, funlabel=c("3-year Desease free survival", "5-year Desease free survival"), 
                        maxscale=100, 
                        fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.3,0.1))
        #plot
        plot(nom)
        #vali caliberation 
        data.for.nomogram.vali<-ValidationData[,step.variable]
        data.for.nomogram.vali<-data.for.nomogram.vali[complete.cases(data.for.nomogram.vali),]
        f.v<-cph(Surv(data.for.nomogram.vali$DFS,data.for.nomogram.vali$DFSstatus)~ predict(f, newdata=data.for.nomogram.vali), x=T, y=T, surv=T, time.inc=60) 
        #get cindex
        rcorr.cens(predict(f, newdata=data.for.nomogram.vali),Surv(data.for.nomogram.vali$DFS,data.for.nomogram.vali$DFSstatus))
        cal.v<-calibrate(f.v, cmethod="KM", method="boot", u=60, m=40, B=1000)
        plot(cal.v)
        #external validation 
        data.for.nomogram.external.vali<-extervaliData[,step.variable]
        data.for.nomogram.external.vali<-data.for.nomogram.external.vali[complete.cases(data.for.nomogram.external.vali),]
        f.e.v<-cph(Surv(data.for.nomogram.external.vali$DFS,data.for.nomogram.external.vali$DFSstatus)~ predict(f, newdata=data.for.nomogram.external.vali), x=T, y=T, surv=T, time.inc=60) 
        cal.e.v<-calibrate(f.e.v, cmethod="KM", method="boot", u=60, m=60, B=1000)
        #get c-index
        rcorr.cens(predict(f, newdata=data.for.nomogram.external.vali),Surv(data.for.nomogram.external.vali$DFS,data.for.nomogram.external.vali$DFSstatus))
        plot(cal.e.v)
        
        #nomogram for OS
        library(rms)
        step.variable<-c("ageStatus","NI","VI","pN_status","riskscoreStatus","OS","OSstatus")
        data.for.nomogram.os<-trainingData[,step.variable]
        data.for.nomogram.os<-data.for.nomogram.os[complete.cases(data.for.nomogram.os),]
        #package data 
        dd2 <- datadist(data.for.nomogram.os)
        options(datadist="dd2")
        #regression data 
        f.os <- cph(Surv(data.for.nomogram.os$OS, data.for.nomogram.os$OSstatus==1) ~ riskscoreStatus+ageStatus+pN_status+NI+VI, x=T, y=T, surv=T, data=data.for.nomogram.os, time.inc=60)
        surv <- Survival(f.os)
        #nomogram
        nom <- nomogram(f, fun=list(function(x) surv(36, x),function(x) surv(60, x)),  
                        lp=F, funlabel=c("3-year Overall surival", "5-year Overall surival"), 
                        maxscale=100, 
                        fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.3,0.1))
        #plot
        plot(nom)
        #vali caliberation 
        data.for.nomogram.vali<-ValidationData[,step.variable]
        data.for.nomogram.vali<-data.for.nomogram.vali[complete.cases(data.for.nomogram.vali),]
        f.v<-cph(Surv(data.for.nomogram.vali$OS,data.for.nomogram.vali$OSstatus)~ predict(f, newdata=data.for.nomogram.vali), x=T, y=T, surv=T, time.inc=60) 
        #get cindex
        rcorr.cens(predict(f, newdata=data.for.nomogram.vali),Surv(data.for.nomogram.vali$OS,data.for.nomogram.vali$OSstatus))
        cal.v<-calibrate(f.v, cmethod="KM", method="boot", u=60, m=35, B=1000)
        plot(cal.v)
        #external validation 
        data.for.nomogram.external.vali<-extervaliData[,step.variable]
        data.for.nomogram.external.vali<-data.for.nomogram.external.vali[complete.cases(data.for.nomogram.external.vali),]
        f.e.v<-cph(Surv(data.for.nomogram.external.vali$OS,data.for.nomogram.external.vali$OSstatus)~ predict(f, newdata=data.for.nomogram.external.vali), x=T, y=T, surv=T, time.inc=60) 
        cal.e.v<-calibrate(f.e.v, cmethod="KM", method="boot", u=60, m=60, B=1000)
        #get c-index
        rcorr.cens(predict(f, newdata=data.for.nomogram.external.vali),Surv(data.for.nomogram.external.vali$OS,data.for.nomogram.external.vali$OSstatus))
        plot(cal.e.v)

#nomogram ROC----
        #cirScore + nomogram and individual factor        
        
        cox.fit.nomogram.DFS <- cph(Surv(data.for.nomogram$DFS, data.for.nomogram$DFSstatus) ~ riskscoreStatus+ageStatus+pN_status+NI+VI , data=trainingData)
        trainingData$nomo.score<-cox.fit.nomogram.DFS$linear.predictors
        ValidationData$nomo.score<--predict(cox.fit.nomogram.DFS,ValidationData)
        extervaliData$nomo.score<--predict(cox.fit.nomogram.DFS,extervaliData)
        interestVariable.ROC<-c("nomo.score","riskscore","Age","pN_status","NI","VI")
        plotmulti_ROC_with_test(trainingData[,interestVariable.ROC],trainingData$DFSstatus)
        plotmulti_ROC_with_test(ValidationData[,interestVariable.ROC],ValidationData$DFSstatus)
        plotmulti_ROC_with_test(extervaliData[,interestVariable.ROC],extervaliData$DFSstatus)
        


