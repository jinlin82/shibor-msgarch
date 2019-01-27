library(zoo)
library(forecast)
library(tseries)
library(xts)
library(ggplot2)
library(FinTS)
library(moments)
library(devtools)
library(MSGARCH)
library(fBasics)
library(timeDate)
library(timeSeries)
library(fGarch)
library(stats)
library(knitr)
rm(list=ls())
path <- "E:/github_repo/shibor-msgarch/data" 
setwd(path)
data=read.csv("Shibor_data.csv",header=T)
shibor<-data[,2]
shibor.date<-as.Date(data[,1])

#计算收益率
shibor.rt<-diff(log(shibor))

#构建残差序列
fit=arima(shibor.rt,order = c(1,0,2))
shibor.rt.r<-stats::residuals(fit,standardize=F)


#输出各regime=1,regime=2,regime=3时的最优模型
spec1=CreateSpec(variance.spec = list(model=c("sGARCH")),
                 distribution.spec = list(distribution=c("norm")))
spec2=CreateSpec(variance.spec = list(model=c("sGARCH","sGARCH")),
                 distribution.spec = list(distribution=c("sstd","ged")))
spec3=CreateSpec(variance.spec = list(model=c("sGARCH","sGARCH","sGARCH")),
                 distribution.spec = list(distribution=c("sged","std","sged")))
models=list(spec1,spec2,spec3)


### VaR 计算

#求VaR
n.ots    <- 700 # number of out-of-sample evaluation 
n.its    <- 2359 # fit sample size
alpha    <- 0.05 # risk Level 
k.update <- 100  # estimation frequency

## Initialization 
VaR   <- matrix(NA, nrow = n.ots, ncol = length(models))
y.ots <- matrix(NA, nrow = n.ots, ncol = 1)
model.fit <- vector(mode = "list", length =length(models))
# iterate over out-of-sample time
for (i in 1:n.ots) {
  cat("Backtest - Iteration: ", i, "\n")
  y.its    <- shibor.rt.r[i:(n.its + i - 1)] # in-sample data 
  y.ots[i] <- shibor.rt.r[n.its + i]         # out-of-sample data
  for (j in 1:length(models)){
    if (k.update == 1 || i %% k.update == 1) {
      cat("Model is reestimated\n")
      model.fit[[j]] <- FitML(spec = models[[j]], data = y.its, 
                              ctr = list(do.se = FALSE)) 
    }
    
    # calculate VaR 1-step ahead
    VaR[i,j] <- Risk(model.fit[[j]]$spec, par = model.fit[[j]]$par,
                     data = y.its,
                     n.ahead = 1,
                     alpha   = alpha,
                     do.es   = FALSE,
                     do.its  = FALSE)$VaR
  }
  
  v=data.frame(VaR,y.ots) 
}
colnames(v)=c("MS(1)-norm","MS(2)-sstd-ged","MS(3)-sged-std-sged","实际值")
kable(v,caption = "VaR预测值与实际值对比表")


library("zoo")
time.index <- zoo::index(shibor.rt.r)[(n.its + 1):(n.ots + n.its)]
y_ots <- zoo::zoo(y.ots, order.by = time.index)
VaR   <- zoo::zoo(VaR, order.by = time.index)

par(mfrow = c(1, 1))
plot(y_ots, type = 'p', las = 1, xlab = "Date (year)",
     ylab = "", col = "black", cex.axis = 1.5, cex.lab = 1.5,cex=0.7, pch = 19)
lines(VaR[,1], type = 'l', col = "red", lwd = 2, lty = "dashed")
lines(VaR[,2], type = 'l', col = "green", lwd = 2, lty = "dashed")
lines(VaR[,3], type = 'l', col = "blue", lwd = 2)
legend("topleft", legend = c("VaR 5% - MS(1)-norm", "VaR 5% - MS(2)-sstd-ged","VaR 5% - MS(3)-sged-std-sged"), 
       col = c("red","green","blue"), lwd = 3, cex = 1, lty = c("dashed", "solid"))
abline(h = 0)#在y=0添加次要刻度线
title("5%显著性水下下的VaR", cex.main = 1.5)

library("GAS")
CC.pval <- DQ.pval <- vector("double", length(models))
for (j in 1:length(models)) {
  test <- GAS::BacktestVaR(data = y.ots, VaR = VaR[,j],
                           alpha = alpha)
  CC.pval[j] <- test$LRcc[2]
  DQ.pval[j] <- test$DQ$pvalue
}

test.print=data.frame(CC.pval,DQ.pval)
#colnames(test.print)=c("MS(1)-norm","MS(2)-sstd-ged","MS(3)-sged-std-sged")
#rownames(test.print)=c("CC.pval","DQ.pval")
kable(test.print,caption = "VaR回测检验")
