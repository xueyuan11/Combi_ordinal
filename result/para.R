library(xlsx)
setwd("C:/Users/o0/Desktop/ordinal data/simulation/Combi_ordinal/result/")
data<-read.xlsx("result_5_4.xlsx",1)
re=data[4:14,3:7]
apply(re,1,function(x){paste(collapse = " & ",x)})
write.table(re, file = "f.txt", sep = "&", row.names=FALSE, col.names=FALSE,quote = F)
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
# alphay = c(-1, 0, 1)
# betay = -.5
# alphax = c(-4:4)
# betax = 1
# eta0 = rep(0,10)
# eta1 = 0.2 * (-4:5)
# eta2= c(-.65,-.54,-.3, .18, .20, .22, .24,.34,.36,.45)
# eta3 = 0.1 * c(-2,0,2,0,-2,-2,0,2,0,-2)
# X_7,Y_4
alphay = c(-1, 0, 1)
betay = -.5
alphax = c(-2:3)
betax = 1
eta0 = rep(0,7)
eta1 = 0.2 * (-3:3)
eta2= c(-.3, .18, .20, .22, .24,.34,.45)
eta3 = 0.1 * c(0,2,0,-2,-2,0,2)
paste(collapse = " , ",eta3)
