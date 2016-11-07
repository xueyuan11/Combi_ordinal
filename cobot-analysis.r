#### Data analysis programms for testing for association between two
#### ordinal variables while adjusting for covariates.
####
#### Chun Li
#### October 20, 2009

#### Chun Li and Bryan E. Shepherd (2009) Test of Association between
#### Two Ordinal Variables while Adjusting for Covariates. JASA (in press)


#### This implementation of our methods assumes proportional odds
#### relationship between y and Z and between x and Z.  In principle,
#### any other multinomial regression analysis can be used in place of
#### proportional odds models.


#### R library requirement:
## Design package is needed for the lrm() function.


#### functions defined:
## diagn():           a modification of diag(), for consistent coding
## tau():             Goodman-Kruskal's gamma
## COBOT.stat():      called by COBOT.emp()
## COBOT.emp():       calculates p-values using empirical distributions
## ordinal.scores():  called by COBOT.scores()
## COBOT.scores():    calculates p-values using asymptotics


#### Data analysis functions:
## COBOT.emp(data, Nemp=1000)
## COBOT.scores(data)
##
## COBOT.emp() calculates p-values using empirical distributions
## COBOT.scores() calculates p-values using asymptotics
##
## The argument "data" must be a list with three elements:
##   y: a vector with integer values representing levels (1 to ny)
##   x: a vector with integer values representing levels (1 to nx)
##   z: a vector or matrix of numerical covariate values
## The lengths of y, x, z (or nrow(z) if matrix) must be the same.
##
## Subjects with missing data are removed by these functions.
## 
## Variables y and x are recoded to so that if variable x (y) has nx (ny)
## categories, the categories are coded from 1 to nx (ny).
##
## For COBOT.emp(), two sets of p-values are estimated:
##   Set a: proportion of abs(emp) >= abs(test.stat)
##   Set b: two times the proportion of smaller tail



#### This works like diag() except when x is a single integer value.
#### When used on the left hand side, use diag().
diagn = function(x) {
  diag(x, length(x), length(x))
}


#### Goodman-Kruskal's gamma.  Unfortunately we had started using "tau"
#### all over the place before realizing this.
tau = function(M) {
  nrow = nrow(M)
  ncol = ncol(M)
  scon = sdis = 0
  for(i in 1:(nrow-1)) {
    for(j in 1:ncol) {
      if(j<ncol)
        scon = scon + M[i,j] * sum(M[(i+1):nrow, (j+1):ncol])
      if(j>1)
        sdis = sdis + M[i,j] * sum(M[(i+1):nrow, 1:(j-1)])
    }
  }

  list(scon=scon, sdis=sdis, tau = (scon-sdis)/(scon+sdis))
}


#### This function is called from COBOT.emp().
#### Calculation of the three test statistics
COBOT.stat = function(data, nx0, ny0) {
  require(rms)  ## for the lrm() function

  ## Ensure code works if called from somewhere else (not COBOT.emp()).
  ## Make z a matrix if it is a vector.  This makes later coding consistent.
  if(!is.matrix(data$z)) data$z = matrix(data$z, ncol=1)

  x = data$x
  y = data$y
  z = data$z

  nx = length(table(x))
  ny = length(table(y))

  nax = nx - 1
  nay = ny - 1
  nz = ncol(z)
  N = length(y)

  ## Fit separate proportional odds models under the null
  ## lrm() may fail if the slope is very close to zero.  This is rare but
  ## needs to be taken care of.  NAs are returned when this happens.
  ## polr() in MASS library seems to have problems too.
##  modxz = lrm(x ~ z)
##  modyz = lrm(y ~ z)
  assign("EE", 0, pos=1)
  suppressWarnings(tryCatch(assign("modxz", lrm(x ~ z, tol=1e-50, maxit=100)),
                            error = function(x) assign("EE", 1, pos=1)))
  suppressWarnings(tryCatch(assign("modyz", lrm(y ~ z, tol=1e-50, maxit=100)),
                            error = function(x) assign("EE", 1, pos=1)))
  if(EE == 1) return(list(stat=rep(NA,1), pred.probx=NULL, pred.proby=NULL))

  ## Obtain predicted cummulative probabilities.
  ## Each column is for a subject
  px = (1 + exp(outer(modxz$coeff[1:nax],
                      as.vector(z%*%modxz$coeff[-(1:nax)]), "+"))) ^ (-1)
  py = (1 + exp(outer(modyz$coeff[1:nay],
                      as.vector(z%*%modyz$coeff[-(1:nay)]), "+"))) ^ (-1)

  ## sum of predicted probability matrix
  pred.probx = rbind(px,1) - rbind(0,px)
  pred.proby = rbind(py,1) - rbind(0,py)
  sum.matrix = tcrossprod(pred.proby, pred.probx)

  ## T1 is to contrast tau of observed and expected frequencies.
  T1 = tau(table(y,x))$tau - tau(sum.matrix)$tau

  ## Fill in the categories that are missed in the empirical data set
  if(nx0 > nx) {
    px.fill = matrix(0, nx0-1, N)
    catx = as.numeric(names(table(x)))
    kk = 0
    for(i in 1:(nx0-1)) {
      if(i %in% catx) kk = kk + 1
      if(kk == 0) {
        px.fill[i,] = 0
      } else if(kk == nx) {
        px.fill[i,] = 1
      } else px.fill[i,] = px[kk,]
    }
    px = px.fill
  }

  if(ny0 > ny) {
    py.fill = matrix(0, ny0-1, N)
    caty = as.numeric(names(table(y)))
    kk = 0
    for(i in 1:(ny0-1)) {
      if(i %in% caty) kk = kk + 1
      if(kk == 0) {
        py.fill[i,] = 0
      } else if(kk == ny) {
        py.fill[i,] = 1
      } else py.fill[i,] = py[kk,]
    }
    py = py.fill
  }

  ## components of residuals
  low.x = rbind(0,px)[cbind(x,1:N)]
  hi.x = 1 - rbind(px,1)[cbind(x,1:N)]
  low.y = rbind(0,py)[cbind(y,1:N)]
  hi.y = 1 - rbind(py,1)[cbind(y,1:N)]

  ## T2 is correlation of residuals
  T2 = cor(hi.x - low.x, hi.y - low.y)

  Cprob = low.x*low.y + hi.x*hi.y
  Dprob = low.x*hi.y + hi.x*low.y
  Cmean = mean(Cprob)
  Dmean = mean(Dprob)

  ## T3 is Cmean-Dmean
  T3 = Cmean - Dmean
  modyzx.con = lrm(y ~ z + x, data=data)
  modyzx.dis = lrm(y ~ z + x,data=list(y=data$y,z=data$z,x=as.factor(data$x)))
  pval.con = anova(modyzx.con)[2,3]
  pval.dis = anova(modyzx.dis)[2,3]
  
  T4=-2*log(pval.con)-2*log(pval.dis)
  w=1/sqrt(nx-1)
  T5=-2*log(pval.con)*(1-w)-2*w*log(pval.dis)
  
  list(stat=c(T1,T2,T3,T4,T5),pval=c(pval.con,pval.dis),
       pred.probx=pred.probx,
       pred.proby=pred.proby)
}


#### P-value estimation for the three test statistics using empirical
#### distributions.  Two sets of p-values are estimated:
####   Set a: proportion of abs(emp) >= abs(test.stat)
####   Set b: two times the proportion of smaller tail
COBOT.emp = function(data, Nemp=1000) {
  ## Make z a matrix if it is a vector.  This makes later coding consistent.
  if(!is.matrix(data$z)) data$z = matrix(data$z, ncol=1)

  ## Check missing data.
  exclude = is.na(data$y) | is.na(data$x) |
            apply(data$z, 1, function(x) sum(is.na(x)))

  ## Exclude subjects with missing data.
  ## Ensure categories are coded from 1 to ny and from 1 to nx.
  data$y = as.numeric(as.factor(data$y[!exclude]))
  data$x = as.numeric(as.factor(data$x[!exclude]))
  data$z = data$z[!exclude, ]

  ## This is ensure data$z is a matrix.  The above statement will make
  ## data$z a vector if ncol=1.
  if(!is.matrix(data$z)) data$z = matrix(data$z, ncol=1)

  nx = length(table(data$x))
  ny = length(table(data$y))
  nxny = nx*ny
  nxnym1 = nxny - 1
  N = length(data$y)

  test.stat = COBOT.stat(data, nx, ny)

  ## cum.mult.prob has cummalative probabilities
  cum.mult.prob = matrix(, nxnym1, N)
  for(i in 1:N)
    cum.mult.prob[,i] = cumsum(outer(test.stat$pred.proby[,i], test.stat$pred.probx[,i]))[-nxny]

  ## Generate empirical distribution of the test statistic under the null
  ## First, generate data using the product predicted probabilities
  xyemp = matrix(,N,Nemp)
  aa = matrix(runif(N*Nemp), N)
  for(k in 1:N)
    xyemp[k,] = outer(aa[k,], cum.mult.prob[,k], ">") %*% rep(1,nxnym1)
  xemp = xyemp %/% ny + 1
  yemp = xyemp %% ny + 1

  test.stat.emp = matrix(,5,Nemp)
  for(j in 1:Nemp) {
    test.stat.emp[,j] = COBOT.stat(list(x=xemp[,j], y=yemp[,j], z=data$z),
                   nx, ny)$stat
  }
#  print(sum(is.na(test.stat.emp[1,])))

  ## Obtain p-values.  Since test.stat.emp may have NAs, apply() is
  ## used instead of the faster matrix multiplication.
#  pvala = ((abs(test.stat$stat) <= abs(test.stat.emp)) %*% rep(1,Nemp))/Nemp
#  tmp = ((test.stat$stat <= test.stat.emp) %*% rep(1,Nemp))/Nemp
  pvala = apply(abs(test.stat$stat) <= abs(test.stat.emp), 1, mean, na.rm=T)
  c(pvala,test.stat$pval)
}


#### This function is called from COBOT.scores().
#### It calculates all values needed for estimating equations.
ordinal.scores = function(y, X) {
  ## y is a numeric vector
  ## X is a vector or matrix with one or more columns.
  require(rms)  ## for the lrm() function

  ## Ensure code works if called from somewhere else (not COBOT.scores()).
  ## Make X a matrix if it is a vector.  This makes later coding consistent.
  if(!is.matrix(X)) X = matrix(X, ncol=1)

  ## N: number of subjects; ny: number of y categories
  N = length(y)
  ny = length(table(y))

  ## na, nb: number of parameters in alpha and beta
  na = ny - 1
  nb = ncol(X)
  npar = na + nb
  ## Z is defined as in McCullagh (1980) JRSSB 42:109-142
  Z = outer(y, 1:ny, "<=")

  ## Fit proportional odds model and obtain the MLEs of parameters.
  mod = lrm(y ~ X, tol=1e-50, maxit=100)
  alpha = -mod$coeff[1:na]
  beta = -mod$coeff[-(1:na)]

  ## Scores are stored for individuals.
  dl.dtheta = matrix(NA, N, npar)
  ## Information matrices are stored as sums over all individuals.
  d2l.dalpha.dalpha = matrix(0,na,na)
  d2l.dalpha.dbeta = matrix(0,na,nb)
  d2l.dbeta.dbeta = matrix(0,nb,nb)
  d2l.dbeta.dalpha = matrix(0,nb,na)
  ## Predicted probabilities p0 and dp0.dtheta are stored for individuals.
  p0 = matrix(,N,ny)
  dp0.dtheta = array(0,c(N,ny,npar))
  ## Cumulative probabilities
  Gamma = matrix(0,N,na)
  dgamma.dtheta = array(0,c(N,na,npar))

  for (i in 1:N) {
    z = Z[i,]  ## z has length ny
    x = X[i,]

    ## gamma and phi are defined as in McCullagh (1980)
    gamma = 1 - 1/(1 + exp(alpha + sum(beta*x))) ## gamma has length na
    diffgamma = diff(c(gamma,1))
    invgamma = 1/gamma
    invgamma2 = invgamma^2
    invdiffgamma = 1/diffgamma
    invdiffgamma2 = invdiffgamma^2
    phi = log(gamma / diffgamma) ## phi has length na
    Gamma[i,] = gamma


#### Some intermediate derivatives
    ## g(phi) = log(1+exp(phi))
    dg.dphi = 1 - 1/(1 + exp(phi))
    ## l is the log likelihood (6.3) in McCullagh (1980)
    dl.dphi = z[-ny] - z[-1] * dg.dphi
    t.dl.dphi = t(dl.dphi)

    ## dphi.dgamma is a na*na matrix with rows indexed by phi
    ## and columns indexed by gamma
    dphi.dgamma = matrix(0,na,na)
    diag(dphi.dgamma) = invgamma + invdiffgamma
    if(na > 1)
      dphi.dgamma[cbind(1:(na-1), 2:na)] = -invdiffgamma[-na]

    dgamma.base = gamma * (1-gamma)
    dgamma.dalpha = diagn(dgamma.base)
    dgamma.dbeta = dgamma.base %o% x
    dgamma.dtheta[i,,] = cbind(dgamma.dalpha, dgamma.dbeta)

    d2gamma.base = gamma * (1-gamma) * (1-2*gamma)

    ##
    d2l.dphi.dphi = diagn(-z[-1] * dg.dphi * (1-dg.dphi))
    d2l.dphi.dalpha = d2l.dphi.dphi %*% dphi.dgamma %*% dgamma.dalpha
    d2l.dphi.dbeta = d2l.dphi.dphi %*% dphi.dgamma %*% dgamma.dbeta

    ##
    d2phi.dgamma.dalpha = array(0,c(na,na,na))
    d2phi.dgamma.dalpha[cbind(1:na,1:na,1:na)] = (-invgamma2 + invdiffgamma2) * dgamma.base
    if(na > 1) {
      d2phi.dgamma.dalpha[cbind(1:(na-1),1:(na-1),2:na)] = -invdiffgamma2[-na] * dgamma.base[-1]
      d2phi.dgamma.dalpha[cbind(1:(na-1),2:na,1:(na-1))] = -invdiffgamma2[-na] * dgamma.base[-na]
      d2phi.dgamma.dalpha[cbind(1:(na-1),2:na,2:na)] = invdiffgamma2[-na] * dgamma.base[-1]
    }

    ##
    d2phi.dgamma.dbeta = array(0,c(na,na,nb))
    rowdiff = matrix(0,na,na)
    diag(rowdiff) = 1
    if(na > 1)
      diag(rowdiff[,-1]) = -1
    d2phi.dgamma.dbeta.comp1 = diagn(-invdiffgamma2) %*% rowdiff %*% dgamma.dbeta
    d2phi.dgamma.dbeta.comp2 = diagn(-invgamma2) %*% dgamma.dbeta - d2phi.dgamma.dbeta.comp1
    for(j in 1:na) {
      d2phi.dgamma.dbeta[j,j,] = d2phi.dgamma.dbeta.comp2[j,]
      if(j < na)
        d2phi.dgamma.dbeta[j,j+1,] = d2phi.dgamma.dbeta.comp1[j,]
    }

    ##
    d2gamma.dalpha.dbeta = array(0,c(na,na,nb))
    for(j in 1:na)
      d2gamma.dalpha.dbeta[j,j,] = d2gamma.base[j] %o% x

    ##
    d2gamma.dbeta.dbeta = d2gamma.base %o% x %o% x


#### First derivatives of log-likelihood (score functions)
    dl.dalpha = dl.dphi %*% dphi.dgamma %*% dgamma.dalpha
    dl.dbeta = dl.dphi %*% dphi.dgamma %*% dgamma.dbeta
    dl.dtheta[i,] = c(dl.dalpha, dl.dbeta)


#### Second derivatives of log-likelihood
#### Since first derivative is a sum of terms each being a*b*c,
#### second derivative is a sum of terms each being (a'*b*c+a*b'*c+a*b*c').

#### d2l.dalpha.dalpha
    ## Obtain aprime.b.c
    ## Transpose first so that matrix multiplication is meaningful.
    ## Then transpose so that column is indexed by second alpha.
    aprime.b.c = t(crossprod(d2l.dphi.dalpha, dphi.dgamma %*% dgamma.dalpha))

    ## Obtain a.bprime.c
    ## run through the index of second alpha
    a.bprime.c = matrix(,na,na)
    for(j in 1:na)
      a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dalpha[,,j] %*% dgamma.dalpha

    ## Obtain a.b.cprime
    ## cprime = d2gamma.dalpha.dalpha = 0 if indices of the two alphas differ.
    d2gamma.dalpha.dalpha = diagn(d2gamma.base)
    a.b.cprime = diagn(as.vector(dl.dphi %*% dphi.dgamma %*% d2gamma.dalpha.dalpha))

    ## summing over individuals
    d2l.dalpha.dalpha = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dalpha.dalpha


#### d2l.dalpha.dbeta
    aprime.b.c = t(crossprod(d2l.dphi.dbeta, dphi.dgamma %*% dgamma.dalpha))
    a.bprime.c = a.b.cprime = matrix(,na,nb)
    for(j in 1:nb) {
      a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dbeta[,,j] %*% dgamma.dalpha
      a.b.cprime[,j] = t.dl.dphi %*% dphi.dgamma %*% d2gamma.dalpha.dbeta[,,j]
    }
    d2l.dalpha.dbeta = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dalpha.dbeta


#### d2l.dbeta.dalpha
#    dl.dbeta = dl.dphi %*% dphi.dgamma %*% dgamma.dbeta
#    aprime.b.c = t(crossprod(d2l.dphi.dalpha, dphi.dgamma %*% dgamma.dbeta))
#    a.bprime.c = a.b.cprime = matrix(,na,nb)
#    for(j in 1:nb) {
#      a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dalpha[,,j] %*% dgamma.dbeta
#      a.b.cprime[,j] = t.dl.dphi %*% dphi.dgamma %*% d2gamma.dbeta.dalpha[,,j]
#    }
#    d2l.dbeta.dalpha = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dbeta.dalpha


#### d2l.dbeta.dbeta
    aprime.b.c = t(crossprod(d2l.dphi.dbeta, dphi.dgamma %*% dgamma.dbeta))
    a.bprime.c = a.b.cprime = matrix(,nb,nb)
    for(j in 1:nb) {
      a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dbeta[,,j] %*% dgamma.dbeta
      a.b.cprime[,j] = t.dl.dphi %*% dphi.dgamma %*% d2gamma.dbeta.dbeta[,,j]
    }
    d2l.dbeta.dbeta = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dbeta.dbeta 


#### Derivatives of predicted probabilities
    p0[i,] = diff(c(0, gamma, 1))

    rowdiff = matrix(0,ny,na)
    diag(rowdiff) = 1
    rowdiff[cbind(2:ny,1:na)] = -1
    dp0.dalpha = rowdiff %*% dgamma.dalpha
    dp0.dbeta = rowdiff %*% dgamma.dbeta

    dp0.dtheta[i,,] = cbind(dp0.dalpha, dp0.dbeta)
  }

#### Final assembly
  d2l.dtheta.dtheta = rbind(
    cbind(d2l.dalpha.dalpha, d2l.dalpha.dbeta),
    cbind(t(d2l.dalpha.dbeta), d2l.dbeta.dbeta))

  ## sandwich variance estimate: ABA', where
  ## A = (-d2l.dtheta.dtheta/N)^(-1)
  ## B = B0/N
  ## One way of coding:
  ## A0 = solve(d2l.dtheta.dtheta)
  ## B0 = t(dl.dtheta) %*% dl.dtheta
  ## var.theta = A0 %*% B0 %*% t(A0)
  ## Suggested coding for better efficiency and accuracy
  SS = solve(d2l.dtheta.dtheta, t(dl.dtheta))
  var.theta = tcrossprod(SS, SS)

  ## The sum of scores should be zero at the MLE.
  ##   apply(dl.dtheta, 2, sum)
  ## Sandwich variance estimate should be similar to information matrix, I,
  ## which is the same as the lrm() output mod$var.
  ##   I = -solve(d2l.dtheta.dtheta)
  ##   print(I)
  ##   print(mod$var)
  ##   print(var.theta)

  list(mod = mod,
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       var.theta = var.theta,
       p0 = p0,
       dp0.dtheta = dp0.dtheta,
       Gamma = Gamma,
       dgamma.dtheta = dgamma.dtheta)
}


#### P-value estimation for the three test statistics using M-estimation.
COBOT.scores = function(data) {
  ## Make z a matrix if it is a vector.  This makes later coding consistent.
  if(!is.matrix(data$z)) data$z = matrix(data$z, ncol=1)

  ## Check missing data.
  exclude = is.na(data$y) | is.na(data$x) |
            apply(data$z, 1, function(x) sum(is.na(x)))

  ## Exclude subjects with missing data.
  ## Ensure categories are coded from 1 to ny and from 1 to nx.
  data$y = as.numeric(as.factor(data$y[!exclude]))
  data$x = as.numeric(as.factor(data$x[!exclude]))
  data$z = data$z[!exclude, ]

  ## This is ensure data$z is a matrix.  The above statement will make
  ## data$z a vector if ncol=1.
  if(!is.matrix(data$z)) data$z = matrix(data$z, ncol=1)

  score.xz = ordinal.scores(data$x, data$z)
  score.yz = ordinal.scores(data$y, data$z)

  npar.xz = dim(score.xz$dl.dtheta)[2]
  npar.yz = dim(score.yz$dl.dtheta)[2]
  xx = data$x; yy = data$y
  nx = length(table(xx))
  ny = length(table(yy))
  N = length(yy)


#### Asymptotics for T3 = mean(Cprob) - mean(Dprob)
  ##  If gamma.x[0]=0 and gamma.x[nx]=1, then
  ##  low.x = gamma.x[x-1], hi.x = (1-gamma.x[x])
  low.x = cbind(0, score.xz$Gamma)[cbind(1:N, xx)]
  low.y = cbind(0, score.yz$Gamma)[cbind(1:N, yy)]
  hi.x = cbind(1-score.xz$Gamma, 0)[cbind(1:N, xx)]
  hi.y = cbind(1-score.yz$Gamma, 0)[cbind(1:N, yy)]
  Cprob = low.x*low.y + hi.x*hi.y
  Dprob = low.x*hi.y + hi.x*low.y
  mean.Cprob = mean(Cprob)
  mean.Dprob = mean(Dprob)
  T3 = mean.Cprob - mean.Dprob

  dgammax.dthetax = score.xz$dgamma.dtheta
  dgammay.dthetay = score.yz$dgamma.dtheta
  dlowx.dthetax = dhix.dthetax = matrix(, npar.xz, N)
  dlowy.dthetay = dhiy.dthetay = matrix(, npar.yz, N)
  for(i in 1:N) {
    if (xx[i] == 1) {
      dlowx.dthetax[,i] <- 0
    } else {
      dlowx.dthetax[,i] <- dgammax.dthetax[i,xx[i]-1,]
    }
    
    if (yy[i] == 1) {
      dlowy.dthetay[,i] <- 0
    } else {
      dlowy.dthetay[,i] <- dgammay.dthetay[i,yy[i]-1,]
    }
    
    if (xx[i] == nx) {
      dhix.dthetax[,i] <- 0
    } else {
      dhix.dthetax[,i] <- -dgammax.dthetax[i,xx[i],]
    }
    
    if (yy[i] == ny) {
      dhiy.dthetay[,i] <- 0
    } else {
      dhiy.dthetay[,i] <- -dgammay.dthetay[i,yy[i],]
    }
    
    ## Old incorrect code:
    #dlowx.dthetax[,i] = ifelse(xx[i] == 1, 0, dgammax.dthetax[i,xx[i]-1,])
    #dlowy.dthetay[,i] = ifelse(yy[i] == 1, 0, dgammay.dthetay[i,yy[i]-1,])
    #dhix.dthetax[,i] = ifelse(xx[i] == nx, 0, -dgammax.dthetax[i,xx[i],])
    #dhiy.dthetay[,i] = ifelse(yy[i] == ny, 0, -dgammay.dthetay[i,yy[i],])
  }
  dCsum.dthetax = dlowx.dthetax %*% low.y + dhix.dthetax %*% hi.y
  dCsum.dthetay = dlowy.dthetay %*% low.x + dhiy.dthetay %*% hi.x
  dDsum.dthetax = dlowx.dthetax %*% hi.y + dhix.dthetax %*% low.y
  dDsum.dthetay = dlowy.dthetay %*% hi.x + dhiy.dthetay %*% low.x

  dT3sum.dtheta = c(dCsum.dthetax-dDsum.dthetax, dCsum.dthetay-dDsum.dthetay)

  ## Estimating equations for (theta, p3)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## p3 is the "true" value of test statistic, and the equation is
  ## p3 - (Ci-Di)
  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta, T3-(Cprob-Dprob))

  ## sandwich variance estimate of var(thetahat)
  Ntheta = npar.xz + npar.yz + 1
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[Ntheta, -Ntheta] = -dT3sum.dtheta
  A[Ntheta, Ntheta] = N

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod(SS, SS)
  varT3 = var.theta[Ntheta, Ntheta]

  pvalT3 = 2 * pnorm(-abs(T3)/sqrt(varT3))


#### Asymptotics for T4 = (mean(Cprob)-mean(Dprob))/(mean(Cprob)+mean(Dprob))
  T4 = (mean.Cprob - mean.Dprob)/(mean.Cprob + mean.Dprob)

  ## Estimating equations for (theta, P4)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## P4 is a vector of (cc, dd, p4).  Their corresponding equations are:
  ## cc - Ci
  ## dd - Di
  ## p4 - (cc-dd)/(cc+dd)
  ## Then p4 is the "true" value of test statistic.
  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta,
    mean.Cprob - Cprob, mean.Dprob - Dprob, 0)

  ## sandwich variance estimate of var(thetahat)
  Ntheta = npar.xz + npar.yz + 3
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[Ntheta-3+(1:3), Ntheta-3+(1:3)] = diag(N, 3)
  A[Ntheta-2, 1:(npar.xz+npar.yz)] = -c(dCsum.dthetax, dCsum.dthetay)
  A[Ntheta-1, 1:(npar.xz+npar.yz)] = -c(dDsum.dthetax, dDsum.dthetay)

  revcpd = 1/(mean.Cprob + mean.Dprob)
  dT4.dcpd = (mean.Cprob-mean.Dprob)*(-revcpd^2)
  A[Ntheta, Ntheta-3+(1:2)] = -N * c(revcpd+dT4.dcpd, -revcpd+dT4.dcpd)

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod(SS, SS)
  varT4 = var.theta[Ntheta, Ntheta]

  pvalT4 = 2 * pnorm(-abs(T4)/sqrt(varT4))


#### Asymptotics for T2 = cor(hi.x - low.x, hi.y - low.y)
  xresid = hi.x - low.x
  yresid = hi.y - low.y
  T2 = cor(xresid, yresid)

  xbyyresid = xresid * yresid
  xresid2 = xresid^2
  yresid2 = yresid^2
  mean.xresid = mean(xresid)
  mean.yresid = mean(yresid)
  mean.xbyyresid = mean(xbyyresid)
  ## T2 also equals numT2 / sqrt(varprod) = numT2 * revsvp
  numT2 = mean.xbyyresid - mean.xresid * mean.yresid
  var.xresid = mean(xresid2) - mean.xresid^2
  var.yresid = mean(yresid2) - mean.yresid^2
  varprod = var.xresid * var.yresid
  revsvp = 1/sqrt(varprod)
  dT2.dvarprod = numT2 * (-0.5) * revsvp^3

  ## Estimating equations for (theta, P5)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## P5 is a vector (ex, ey, crossxy, ex2, ey2, p5).
  ## Their corresponding equations are:
  ## ex - (hi.x-low.x)[i]
  ## ey - (hi.y-low.y)[i]
  ## crossxy - ((hi.x-low.x)*(hi.y-low.y))[i]
  ## ex2 - (hi.x-low.x)[i]^2
  ## ey2 - (hi.y-low.y)[i]^2
  ## p5 - (crossxy-ex*ey)/sqrt((ex2-ex^2)*(ey2-ey^2))
  ## Then p5 is the "true" value of test statistic
  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta,
    mean.xresid - xresid, mean.yresid - yresid, mean.xbyyresid - xbyyresid,
    mean(xresid2) - xresid2, mean(yresid2) - yresid2, 0)

  ## sandwich variance estimate of var(thetahat)
  Ntheta = npar.xz + npar.yz + 6
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[Ntheta-6+(1:6), Ntheta-6+(1:6)] = diag(N, 6)

  dxresid.dthetax = dhix.dthetax - dlowx.dthetax
  dyresid.dthetay = dhiy.dthetay - dlowy.dthetay
  bigpartial = rbind(c(dxresid.dthetax %*% rep(1, N), rep(0, npar.yz)),
    c(rep(0, npar.xz), dyresid.dthetay %*% rep(1, N)),
    c(dxresid.dthetax %*% yresid, dyresid.dthetay %*% xresid),
    c(dxresid.dthetax %*% (2*xresid), rep(0, npar.yz)),
    c(rep(0, npar.xz), dyresid.dthetay %*% (2*yresid)))
  A[Ntheta-6+(1:5), 1:(npar.xz+npar.yz)] = -bigpartial

  smallpartial = N *
    c(-mean.yresid * revsvp + dT2.dvarprod * (-2*mean.xresid*var.yresid),
      -mean.xresid * revsvp + dT2.dvarprod * (-2*mean.yresid*var.xresid),
      revsvp,
      dT2.dvarprod * var.yresid,
      dT2.dvarprod * var.xresid)
  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod(SS, SS)
  varT2 = var.theta[Ntheta, Ntheta]

  pvalT2 = 2 * pnorm(-abs(T2)/sqrt(varT2))


#### Asymptotics for T1 = tau - tau0
  ## dtau0/dtheta
  ## P0 is the sum of product predicted probability matrix with dim(nx,ny)
  P0 = crossprod(score.xz$p0, score.yz$p0) / N

  cdtau0 = tau(P0)
  C0 = cdtau0$scon
  D0 = cdtau0$sdis
  ## C0 = sum_{l>j,m>k} {P0[j,k] * P0[l,m]}
  ## D0 = sum_{l>j,m<k} {P0[j,k] * P0[l,m]}
  dC0.dP0 = matrix(,nx,ny)
  dD0.dP0 = matrix(,nx,ny)
  for(i in 1:nx)
    for(j in 1:ny) {
      dC0.dP0[i,j] = ifelse(i>1 & j>1, sum(P0[1:(i-1), 1:(j-1)]), 0) +
        ifelse(i<nx & j<ny, sum(P0[(i+1):nx, (j+1):ny]), 0)
      dD0.dP0[i,j] = ifelse(i>1 & j<ny, sum(P0[1:(i-1), (j+1):ny]), 0) +
        ifelse(i<nx & j>1, sum(P0[(i+1):nx, 1:(j-1)]), 0)
    }

  ## tau0 = (C0-D0)/(C0+D0)
  dtau0.dC0 = 2*D0/(C0+D0)^2
  dtau0.dD0 =-2*C0/(C0+D0)^2

  ## P0 is already a matrix
  dP0.dtheta.x = array(0, c(nx, ny, npar.xz))
  for(j in 1:ny) {
    aa = matrix(0, nx, npar.xz)
    for(i in 1:N)
      aa = aa + score.xz$dp0.dtheta[i,,] * score.yz$p0[i,j]
    dP0.dtheta.x[,j,] = aa/N
    ## simpler but mind-tickling version
    #dP0.dtheta.x[,j,] = (score.yz$p0[,j] %*% matrix(score.xz$dp0.dtheta,N))/N
  }

  dP0.dtheta.y = array(0, c(nx, ny, npar.yz))
  for(j in 1:nx) {
    aa = matrix(0, ny, npar.yz)
    for(i in 1:N)
      aa = aa + score.yz$dp0.dtheta[i,,] * score.xz$p0[i,j]
    dP0.dtheta.y[j,,] = aa/N
  }

  ## dC0.dtheta and dD0.dtheta
  dC0.dtheta.x = as.numeric(dC0.dP0) %*% matrix(dP0.dtheta.x, nx*ny)
  dD0.dtheta.x = as.numeric(dD0.dP0) %*% matrix(dP0.dtheta.x, nx*ny)
  dC0.dtheta.y = as.numeric(dC0.dP0) %*% matrix(dP0.dtheta.y, nx*ny)
  dD0.dtheta.y = as.numeric(dD0.dP0) %*% matrix(dP0.dtheta.y, nx*ny)

  ## dtau0/dtheta
  dtau0.dtheta.x = dtau0.dC0 * dC0.dtheta.x + dtau0.dD0 * dD0.dtheta.x
  dtau0.dtheta.y = dtau0.dC0 * dC0.dtheta.y + dtau0.dD0 * dD0.dtheta.y


  ## dtau/dPa
  ## tau = (C-D)/(C+D)
  Pa = table(data$x, data$y) / N
  cdtau = tau(Pa)
  C = cdtau$scon
  D = cdtau$sdis
  dtau.dC = 2*D/(C+D)^2
  dtau.dD =-2*C/(C+D)^2

  ## Pa[nx,ny] is not a parameter, but = 1 - all other Pa parameters.
  ## Thus, d.Pa[nx,ny]/d.Pa[i,j] = -1.
  ## Also, d.sum(Pa[-nx,-ny]).dPa[i,j] = 1 when i<nx and j<ny, and 0 otherwise.
  ##
  ## In C = sum_{l>j,m>k} {Pa[j,k] * Pa[l,m]}, Pa[i,j] appears in
  ##   Pa[i,j] * XX (minus Pa[nx,ny] if i<nx & j<ny), and in
  ##   sum(Pa[-nx,-ny]) * Pa[nx,ny].
  ## So, dC.dPa[i,j] = XX (minus Pa[nx,ny] if i<nx & j<ny)
  ##                  + d.sum(Pa[-nx,-ny]).dPa[i,j] * Pa[nx,ny]
  ##                  - sum(Pa[-nx,-ny])
  ##                 = XX (with Pa[nx,ny] if present) - sum(Pa[-nx,-ny])
  ##
  ## D = sum_{l>j,m<k} {Pa[j,k] * Pa[l,m]} doesn't contain Pa[nx,ny]
  dC.dPa = matrix(,nx,ny)
  dD.dPa = matrix(,nx,ny)
  for(i in 1:nx)
    for(j in 1:ny) {
      dC.dPa[i,j] = ifelse(i>1 & j>1, sum(Pa[1:(i-1), 1:(j-1)]), 0) +
        ifelse(i<nx & j<ny, sum(Pa[(i+1):nx, (j+1):ny]), 0) - sum(Pa[-nx,-ny])
      dD.dPa[i,j] = ifelse(i>1 & j<ny, sum(Pa[1:(i-1), (j+1):ny]), 0) +
        ifelse(i<nx & j>1, sum(Pa[(i+1):nx, 1:(j-1)]), 0)
    }

  dtau.dPa = dtau.dC * dC.dPa + dtau.dD * dD.dPa
  dtau.dPa = dtau.dPa[-length(dtau.dPa)]  ## remove the last value


  ## Estimating equations for (theta, phi)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## phi is (p_ij) for (X,Y), and the equations are
  ## I{subject in cell (ij)} - p_ij
  phi.Pa = matrix(0, N, nx*ny)
  phi.Pa[cbind(1:N, xx+(yy-1)*nx)] = 1
  phi.Pa = phi.Pa - rep(1,N) %o% as.numeric(Pa)
  phi.Pa = phi.Pa[,-(nx*ny)]

  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta, phi.Pa)

  ## sandwich variance estimate of var(thetahat, phihat)
  Ntheta = npar.xz + npar.yz + nx*ny-1
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[-(1:(npar.xz+npar.yz)), -(1:(npar.xz+npar.yz))] = -diag(N, nx*ny-1)

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  ##SS = solve(A, t(bigphi))
  ##var.theta = SS %*% t(SS)
  ## Or better yet, no need to explicitly obtain var.theta.  See below.

  ## Test statistic T1 = tau - tau0
  T1 = cdtau$tau - cdtau0$tau
  ## dT.dtheta has length nx + ny + nx*ny-1
  dT1.dtheta = c(-dtau0.dtheta.x, -dtau0.dtheta.y, dtau.dPa)

  ## variance of T, using delta method
  ##varT = t(dT.dtheta) %*% var.theta %*% dT.dtheta
  SS = crossprod(dT1.dtheta, solve(A, t(bigphi)))
  varT1 = sum(SS^2)
  pvalT1 = 2 * pnorm(-abs(T1)/sqrt(varT1))

  c(pvalT1, pvalT2, pvalT3)
}
