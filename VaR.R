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
path <- "F:/github_repo/shibor-msgarch/data" 
setwd(path)
data=read.csv("Shibor_data.csv",header=T)
shibor<-data[,2]
shibor.date<-as.Date(data[,1])

#计算收益率
shibor.rt<-diff(log(shibor))

#构建残差序列
fit=arima(shibor.rt,order = c(1,0,2))
shibor.rt.r<-stats::residuals(fit,standardize=F)


#输出最优模型第216个的结果
spec=CreateSpec(variance.spec = list(model=c("sGARCH","sGARCH","sGARCH")),distribution.spec = 
                  list(distribution=c("sged","std","sged")))
m3.217=FitML(spec,shibor.rt.r)
summary(m3.217)

m3
#求VaR
n.ots    <- 700 # number of out-of-sample evaluation 
n.its    <- 2359 # fit sample size
alpha    <- 0.05 # risk Level 
k.update <- 100  # estimation frequency

## Initialization 
VaR   <- matrix(NA, nrow = n.ots, ncol = 1)
y.ots <- matrix(NA, nrow = n.ots, ncol = 1)
model.fit <- vector(mode = "list", length = 1)

# iterate over out-of-sample time
for (i in 1:n.ots) { 
  cat("Backtest - Iteration: ", i, "\n")
  y.its    <- shibor.rt.r[i:(n.its + i - 1)] # in-sample data 
  y.ots[i] <- shibor.rt.r[n.its + i]         # out-of-sample data

    if (k.update == 1 || i %% k.update == 1) {
      cat("Model is reestimated\n")
      model.fit[[1]] <- FitML(spec = spec, data = y.its, 
                              ctr = list(do.se = FALSE)) 
    }
    # calculate VaR 1-step ahead
    VaR[i,1] <- Risk(model.fit[[1]]$spec, par = model.fit[[1]]$par,
                     data = y.its,
                     n.ahead = 1,
                     alpha   = alpha,
                     do.es   = FALSE,
                     do.its  = FALSE)$VaR
  }                                


## Test the VaR
#install.packages("GAS")
library("GAS")
#for (j in 1:length(models)) { 
test <- GAS::BacktestVaR(data  = y.ots,
                         VaR   = VaR[,1],
                          alpha = alpha)

  CC.pval<- test$LRuc[2]      
  DQ.pval<- test$DQ$pvalue  
#}
#names(CC.pval) <- names(DQ.pval) <- c("GJR-std", "MS2-GJR-std")
?BacktestVaR
print(CC.pval)
print(DQ.pval)



library("zoo")
time.index <- zoo::index(shibor.rt.r)[(n.its + 1):(n.ots + n.its)]
y_ots <- zoo::zoo(y.ots, order.by = time.index)
VaR   <- zoo::zoo(VaR, order.by = time.index)

#pdf(file = "figure6.pdf", height = 13, width = 16, compress = TRUE)
par(mfrow = c(1, 1))
plot(y_ots, type = 'p', las = 1, cex = 0.7, xlab = "Date (year)",ylim=c(-0.4,0.4),
     ylab = "", col = "black", cex.axis = 1.5, cex.lab = 1.5, pch = 19)
lines(VaR[,1], type = 'l', col = "red", lwd =2, lty = "dashed")

abline(h = 0)
title("Backtesing VaR at 5% risk level", cex.main = 1.5)

k=0
for(i in 1:700){
  if(VaR[i,1]>=shibor.rt.r[(2359+i)]){
    k=k+1
  }
}
k
k/700
