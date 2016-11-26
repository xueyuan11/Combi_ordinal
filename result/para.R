library(xlsx)
setwd("C:/Users/o0/Desktop/ordinal data/simulation/Combi_ordinal/result/")
  data<-read.xlsx("new/result54-3.xlsx",1)
re=data[3:13,3:7]
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

###X_3Y_4
alphay = c(-1, 0, 1)
betay = -.5
alphax = c(0,1)
betax = 1
eta0 = rep(0,3)
eta1 = 0.2 * (-1:1)
eta2= c(-.3, .18, .20)
eta3 = 0.1 * c(-2,0,2)
paste(collapse = " , ",eta3)

setwd("J:/Onedirve/OneDrive/Documents/result/cob2/X_5_Y_4")
alpha<- read.xlsx(file = "result_alpha.xlsx",1,header = F)
beta1<-read.xlsx(file = "result_b1.xlsx",1,header = F)
beta2<-read.xlsx(file = "result_b2.xlsx",1,header = F)
beta3<-read.xlsx(file = "result_b3.xlsx",1,header = F)

result<-cbind(alpha,beta1,beta2,beta3)


rownames(result)=c("T1emp", "T2emp", "T3emp","Cob2","Cob4","min2","min4",
                   "Y~X linear","Y~X catego","X~Y linear","X~Y catego")
colnames(result)=c("Null","Linear","Nonlinear","Nonmontonic")
write.xlsx(x = result, file = "C:/Users/o0/Desktop/ordinal data/simulation/Combi_ordinal/result/result.xlsx",
           sheetName = "result", row.names =TRUE,col.names = TRUE)
