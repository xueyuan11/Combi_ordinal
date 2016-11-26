
setwd("C:/Users/o0/Desktop/ordinal data/simulation/Combi_ordinal/testspeed")
#### Need all functions defined in COBOT-analysis.r
source("ca.r")



#### Function for simulating data
generate.data = function(k,alphax, betax, alphay, betay, eta, N) {
  z = rnorm(N,0,1)
  x = y = numeric(N)
  
  ## px is an N x length(alphax) matrix.
  ## Each row has the TRUE cummulative probabilities for each subject.
  px = (1 + exp(- outer(alphax, betax*z, "+"))) ^ (-1)
  aa = runif(N)
  for(i in 1:N)
    x[i] = sum(aa[i] > px[,i])
  x = as.numeric(as.factor(x))
  
  
  py = (1 + exp(- outer(alphay, betay*z+eta[x], "+"))) ^ (-1)
  aa = runif(N)
  for(i in 1:N)
    y[i] = sum(aa[i] > py[,i])
  y = as.numeric(as.factor(y))
  ## y = y+1 may have category gaps if there are small probability categories.
  return(c(x,y,z))
}




#### Start simulations
ptm=proc.time()

N = 500
Nemp = 1000
NREPL = 1000
# 1. X_5. Y_4
alphay = c(-1, 0, 1)
betay = -0.5
alphax = c(-1, 0, 1, 2)
betax = 1
eta0 = rep(0,5)
eta1 =c(-0.3,-0.15,0,0.15,0.3)
eta2= c(-0.21,-0.18,0.10,0.22,0.24)
eta3 = c(-0.35,0.4,0.2,0,-0.15)

# #null distribution
data=as.data.frame(sapply(1:NREPL,generate.data,alphax, betax, alphay, betay, eta0, N))
pval=lapply(data,COBOT.emp,Nemp)

proc.time()-ptm


##Computer Type I error and power
library(xlsx)

p=0.05
pvalm<-matrix(unlist(pval),nrow = NREPL,byrow = T)
alpha<-apply(pvalm<p,2,mean,rm.na=T)

# setwd("J:/Onedirve/OneDrive/Documents/result/cob2/X_5_Y_4")
# write.xlsx(x = alpha, file = "result_alpha.xlsx",
#            sheetName = "result", row.names =F,col.names = F)
