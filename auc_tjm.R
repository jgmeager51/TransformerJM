library(survival)
library(tdROC)
library(JMbayes2)
library(nlme)
surv_pred <- read.csv("C:/Users/jgmea/research/transf/TransformerJM/surv_pred_1.csv",header=FALSE)
testdata <- read.csv("C:/Users/jgmea/research/transf/TransformerJM/test_data1.csv")
traindata <- read.csv("C:/Users/jgmea/research/transf/TransformerJM/train_data1.csv")
tmpdata <- read.csv("C:/Users/jgmea/research/transf/TransformerJM/tmp_data1.csv")

temp.time <- read.csv("C:/Users/jgmea/research/transf/TransformerJM/time_tmp1.csv",header=FALSE)[-1,]
temp.event <- read.csv("C:/Users/jgmea/research/transf/TransformerJM/event_tmp1.csv",header=FALSE)[-1,]
train.time <- read.csv("C:/Users/jgmea/research/transf/TransformerJM/time_train1.csv",header=FALSE)
event.time <- read.csv("C:/Users/jgmea/research/transf/TransformerJM/event_train1.csv",header=FALSE)

#test <- AUC(1-surv_pred, temp.event, temp.time, predtimes)
X <- 1-surv_pred
roc<- tdROC(X = X[,1], Y = temp.time, delta = temp.event, tau = 5, span = 0.05,
            nboot = 0, alpha = 0.05, n.grid = 1000, cut.off = 0.5)
roc$AUC

?tdROC
#surv_pred, e_tmp.numpy().astype(int), t_tmp.numpy(),
#                   e_train.numpy(), t_train.numpy(), LT, np.array(pred_windows)



Brier = function(surv, event, time, event_train, time_train, LT, DeltaT){
  #estimate km curve for BS calculation
  train.surv = cbind.data.frame("event"=event_train, "time"=time_train)
  km = survfit(Surv(time, event)~1, data=train.surv)
  survest = stepfun(km$time, c(1, km$surv))
  
  BS = rep(NA, length(DeltaT))
  for(i in 1:length(DeltaT)){
    pt = LT + DeltaT[i]
    N_vali = length(event)
    
    #BRIER SCORE
    D = rep(0, N_vali) 
    D[time<=pt & event==1] = 1
    
    pi = 1-surv[,i]
    
    km_pts = survest(time)/survest(LT)
    W2 <- D/km_pts
    W1 <- as.numeric(time>pt)/(survest(pt)/survest(LT))
    W <- W1 + W2
    
    BS_pts <- W * (D - pi)^2
    BS[i] = sum(na.omit(BS_pts)) / N_vali
  }
  return(BS)
}
library(pracma)

get_integrated <- function(x, times) {
  integrated_value <- trapz(times, x) / (max(times) - min(times))
  return(integrated_value)
}

Brier(surv_pred, as.integer(temp.event), temp.time, event.time, train.time, 1, .5)
