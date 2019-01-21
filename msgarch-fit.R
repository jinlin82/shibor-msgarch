rm(list = ls())
options(digits = 9)
library(zoo)
library(forecast)
library(tseries)
library(xts)
library(ggplot2)
library(FinTS)
library(moments)
library(devtools)
library(FinTS)
library(MSGARCH)
library(fBasics)
library(timeDate)
library(timeSeries)
library(fBasics)
library(fGarch)
library(forecast)
library(stats)


## path <- "F:/github_repo/shibor-msgarch/" 
## setwd(path)
dat=read.csv("Shibor_data.csv",header=T)
shibor<-dat[,2]
shibor.date<-as.Date(dat[,1])

#计算收益率
shibor.rt<-diff(log(shibor))

#构建残差序列
fit=arima(shibor.rt,order = c(1,0,2))
shibor.rt.r<-residuals(fit,standardize=T)


# 生成 loglik AIC BIC 表 函数
msgarchres <- function(dists, data){
### dists 为 分布名称构成的数据框
    n <- nrow(dists)
    k <- ncol(dists)
    res <- list()
    for(i in 1:n){
        spec=CreateSpec(variance.spec = list(model=rep("sGARCH",k)),
                        distribution.spec =
                            list(distribution=as.vector(as.matrix(dists[i,]))))
        try(temp <- FitML(spec,data), silent = T)
        if(exists("temp")){
            res[[i]] <- c(paste("MSGARCH", k, sep = "-"),
                          paste(dists[i,], collapse = "-"),
                          temp$loglik,
                          summary(temp)$AIC,
                          summary(temp)$BIC)
        } else{
            res[[i]] <- c(paste("MSGARCH", k, sep = "-"),
                          paste(dists[i,], collapse = "-"),
                          rep(NA, 3))  
        }
        print(i)
    }

    res <- as.data.frame(do.call(rbind, res))
    res[,3:5] <- apply(res[3:5], 2, as.numeric)
    names(res) <- c("模型", "分布", "Loglik", "AIC", "BIC")
    return(res)
}


distname=c("norm","snorm","std","sstd","ged","sged")
#各种组合方式
one=expand.grid(distname, stringsAsFactors=F)
two=expand.grid(distname,distname,stringsAsFactors=F)
three=expand.grid(distname,distname,distname,stringsAsFactors=F)
four=expand.grid(distname,distname,distname,distname,stringsAsFactors=F)

#区制为1的MS-GARCH
system.time( stage1 <- msgarchres(one, data=shibor.rt.r) )

write.csv(stage1, './result/stage1.csv')

#区制为2的MS-GARCH
system.time( stage2 <- msgarchres(two, data=shibor.rt.r) )

write.csv(stage2, './result/stage2.csv')

## Donot Run
## 区制为3的MS-GARCH
system.time( stage3 <- msgarchres(three, data=shibor.rt.r) )
write.csv(stage3, './result/stage3.csv')

## Donot Run
#区制为4的MS-GARCH
system.time( stage4 <- msgarchres(four, data=shibor.rt.r) )
write.csv(stage4, './result/stage4.csv')

