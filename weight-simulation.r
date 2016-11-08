#### Simulation programms for testing for association between two
#### ordinal variables while adjusting for covariates.
####
#### Chun Li
#### October 15, 2009

#### Chun Li and Bryan E. Shepherd (2009) Test of Association between
#### Two Ordinal Variables while Adjusting for Covariates. JASA (in press)


#### Here we assume proportional odds relationship between y and Z and
#### between x and Z.  In principle, any other multinomial regression
#### analysis can be used in place of proportional odds models.


#### R library requirement:
## rms package is needed for the lrm() function.


#### functions defined in addition to those in COBOT-analysis.r :
## isotonic():        proportional odds model with an isotonic procedure
## generate.data():   data generation
## ordinalsim():      main function of simulation


#### Data requirements: [COBOT.stat(), COBOT.emp(), COBOT.scores(), isotonic()]
## The argument "data" must be a list with three elements:
##   y: a vector with integer values representing levels (1 to ny)
##   x: a vector with integer values representing levels (1 to nx)
##   z: a vector or matrix of numerical covariate values
## The lengths of y, x, z (or nrow(z) if matrix) must be the same.
#S1 
setwd("C:/Users/o0/Desktop/ordinal data/simulation/Combi_ordinal")
#### Need all functions defined in COBOT-analysis.r
source("cobot-analysis.r")

ptm<-proc.time()
#### Isotonic regression
isotonic = function(data) {
  require(rms)

  mod = lrm(y ~ z + x,data=list(y=data$y,z=data$z,x=as.factor(data$x)))

  xdata = as.factor(data$x)
  xlevels = length(levels(xdata))
  coeff = c(0, rev(rev(mod$coeff)[1:(xlevels-1)]))
  xcoeff = coeff[as.numeric(xdata)]

  #### R's isoreg() function only assumes increasing
  #### The iKnots output of the function can give incorrect results
  iso = iso1 = isoreg(xdata, xcoeff)
  iso2 = isoreg(xdata, -xcoeff)
  if(length(unique(iso1$yf)) < length(unique(iso2$yf))) iso = iso2

  #### Define new x categories: if no jump, all into a single category
  jumppts = which(diff(iso$yf) > 0)
  newx = rep(0, length(xdata))
  if(length(jumppts) > 0) {
    xrank = rank(xdata, ties.method="min")
    for(i in 1:length(jumppts))
      newx = newx + (xrank > jumppts[i])
  }

  newmod = lrm(y ~ z + x,data=list(y=data$y,z=data$z,x=as.factor(newx)))
  anova(newmod)[2,3]
}


#### Function for simulating data
generate.data = function(alphax, betax, alphay, betay, eta, N) {
  z = rnorm(N,0,1)
  x = y = numeric(N)
  
  ## px is an N x length(alphax) matrix.
  ## Each row has the TRUE cummulative probabilities for each subject.
  px = (1 + exp(- outer(alphax, betax*z, "+"))) ^ (-1)
  aa = runif(N)
  for(i in 1:N)
    x[i] = sum(aa[i] > px[,i])
  x = as.numeric(as.factor(x))
  ## x = x+1 may have category gaps if there are small probability categories.
  
  py = (1 + exp(- outer(alphay, betay*z+eta[x], "+"))) ^ (-1)
  aa = runif(N)
  for(i in 1:N)
    y[i] = sum(aa[i] > py[,i])
  y = as.numeric(as.factor(y))
  ## y = y+1 may have category gaps if there are small probability categories.
  
  return(list(x=x, y=y, z=z))
}
  

#### Function for simulation
ordinalsim = function(alphax, betax, alphay, betay, eta, N, Nemp, NREPL) {
  require(rms)

  pval.emp = matrix(,NREPL,7)
  pval.score = matrix(,NREPL,3)
   pval.spl = pval.iso = NULL

  for (ii in 1:NREPL) {
    data = generate.data(alphax, betax, alphay, betay, eta, N)
    pval.emp[ii,] = COBOT.emp(data, Nemp)
    pval.score[ii,] = COBOT.scores(data)
    modspline = lrm(y ~ z + rcs(x, 3), data=data)
    pval.spl[ii] = anova(modspline)[2,3]
   pval.iso[ii] = isotonic(data)
  }

   
   param = list(alpha.yz=alphay, beta.yz=betay,
    alpha.xz=alphax, beta.xz=betax, eta=eta,
    N=N, Nemp=Nemp, Nrepl=NREPL)
  list(par= param, pval = data.frame(pval.emp, pval.score, pval.spl, pval.iso))
}


#### Start simulations


N = 50
Nemp = 100
NREPL = 10
# # 1. X_5. Y_4
# alphay = c(-1, 0, 1)
# betay = -.5
# alphax = c(-1, 0, 1, 2)
# betax = 1
# eta0 = rep(0,5)
# eta1 = 0.2 * (-2:2)
# eta2= c(-.3, .18, .20, .22, .24)
# eta3 = 0.1 * c(-2,0,2,0,-2)
#2. X_10. Y_4
alphay = c(-1, 0, 1)
betay = -.5
alphax = c(-4:4)
betax = 1
eta0 = rep(0,10)
eta1 = 0.2 * (-4:5)
eta2= c(-.65,-.54,-.3, .18, .20, .22, .24,.34,.36,.45)
eta3 = 0.1 * c(-2,0,2,0,-2,-2,0,2,0,-2)

# #null distribution
sim0s = ordinalsim(alphax, betax, alphay, betay, eta, N, Nemp, NREPL)


#linear

sim1s = ordinalsim(alphax, betax, alphay, betay, eta2, N, Nemp, NREPL)

#nonlinear and monoto

sim2s = ordinalsim(alphax, betax, alphay, betay, eta3, N, Nemp, NREPL)

#non-monoto

sim3s = ordinalsim(alphax, betax, alphay, betay, eta, N, Nemp, NREPL)




##Computer Type I error and power
library(xlsx)

p=0.05

alpha<-apply(sim0s$pval<p,2,mean,rm.na=T)
beta=matrix(0,12,1)
beta<-apply(sim1s$pval<p,2,mean,rm.na=T)
beta[,2]<-apply(sim2s$pval<p,2,mean,rm.na=T)
beta[,3]<-apply(sim3s$pval<p,2,mean,rm.na=T)

result<-cbind(alpha,beta)


rownames(alpha)=c("T1emp", "T2emp", "T3emp","X linear","X catego","CobT1","CobT2", "T1s", "T2s","T3s","iso","Spline")
colnames(result)=c("Null","Linear","Nonlinear","Nonmontonic")

setwd("J:/Onedirve/OneDrive/Documents/result/cob/X_10_Y_4")
write.xlsx(x = alpha, file = "result_10_4.xlsx",
           sheetName = "X_10_Y_4", row.names =TRUE,col.names = TRUE)

# apply(result,1,function(x){paste(collapse = " & ",x)})
# write.table(result, file = "foo.txt", sep = "&", row.names=FALSE, col.names=FALSE,
