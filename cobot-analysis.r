

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

#### This works like diag() except when x is a single integer value.
#### When used on the left hand side, use diag().


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

  # ## Fill in the categories that are missed in the empirical data set
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
  
  ##
  modyzx.con = lrm(y ~ z + x, data=data)
  modyzx.dis = lrm(y ~ z + x,data=list(y=data$y,z=data$z,x=as.factor(data$x)))
  pval.conyx = anova(modyzx.con)[2,3]
  pval.disyx = anova(modyzx.dis)[2,3]
  ## conbime 
  T4=-2*log(pval.conyx)-2*log(pval.disyx)
  T5=min(pval.disyx,pval.conyx)
  ##
  modxzy.con = lrm(x ~ z + y, data=data)
  modxzy.dis = lrm(x ~ z + y,data=list(x=data$x,z=data$z,y=as.factor(data$y)))
  pval.conxy = anova(modxzy.con)[2,3]
  pval.disxy = anova(modxzy.dis)[2,3]
  T6=T4-2*log(pval.conxy)-2*log(pval.disxy)
  T7=min(T5,pval.disxy,pval.conxy)
  
  ##
  list(stat=c(T1,T4,T5,T6,T7),pval=c(pval.conyx,pval.disyx,pval.conxy,pval.disxy),
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

  pvala = apply(abs(test.stat$stat) <= abs(test.stat.emp), 1, mean, na.rm=T)
  c(pvala,test.stat$pval)
}

