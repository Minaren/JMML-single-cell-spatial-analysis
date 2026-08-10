rm(list = ls())
mydata <- read.table('mydata.txt',header = T,row.names = 1,sep = '\t')
mydata <- apply(mydata, 2, function(x){as.numeric(x)})

kuan2chang<- function(mydata){
  myresult <- data.frame()
  mytemp <- c()
  for (i in 1:ncol(mydata)) {
    if(i<=ncol(mydata)/2){
      fre <- mydata[,i]
      sample <- rep(colnames(mydata[i]),length(mydata[,i]))
      mic <- rownames(mydata)
      group <- rep('group1',length(mydata[,i]))
      mytemp<-rbind(fre,sample,group,mic)
      mytemp <- t(mytemp)
    }else {
      fre <- mydata[,i]
      sample <- rep(colnames(mydata[i]),length(mydata[,i]))
      group <- rep('group2',length(mydata[,i]))
      mic <- rownames(mydata)
      mytemp<-rbind(fre,sample,group,mic)
      mytemp <- t(mytemp)
    }
    myresult <- rbind(myresult,mytemp)
  }
  return(myresult)
}#定义数据格式转换函数
myresult <- kuan2chang(mydata)#数据格式转换函

biomambatest<-function(mydata){
  mytest <- shapiro.test(mydata)
  if(mytest$p.value >= 0.05){print('样本符合正态分布，请做方差检验')
    mybartest <- bartlett.test(fre~group,data=myresult)
    if(mybartest$p.value>=0.05){print("方差齐性，请做参数检验（t.test）")}else{print('方差非齐性，请做ANOVA检验')}
  }else{print('样本不符合正态分布，请使用krukal.test或wilcoxon')}
}#定义判断函数

biomambatest(mydata)#检验

set.seed(90)
x <- rbinom(15,8,0.7)#生成一个符合正态分布的数据集
biomambatest(x)