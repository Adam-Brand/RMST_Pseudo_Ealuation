######## Analysis of results


sum.func <- function(result, scen){
     
          res.mean <- colMeans(result, na.rm=TRUE)
          
          reject.psu.sand <- sum(result[,5]<0.05 & !is.na(result[,5]))/sum(!is.na(result[,5]))
          reject.psu.ajk <- sum(result[,10]<0.05 & !is.na(result[,10]))/sum(!is.na(result[,10]))
          reject.km <- sum(result[,15]<0.05 & !is.na(result[,15]))/sum(!is.na(result[,15]))
          reject.rp.nai <- sum(result[,20]<0.05 & !is.na(result[,20]))/sum(!is.na(result[,20]))
          reject.rp.asy <- sum(result[,25]<0.05 & !is.na(result[,25]))/sum(!is.na(result[,25]))
          if (ncol(result)==35){
               reject.cox <- sum(result[,30]<0.05 & !is.na(result[,30]))/sum(!is.na(result[,30]))
               reject.ph <- sum(result[,33]<0.05 & !is.na(result[,33]))/sum(!is.na(result[,33]))
               bias.psu.sand <- result[,1]-result[,35]
               bias.psu.ajk <- result[,6]-result[,35]
               bias.km <- result[,11]-result[,35]
               bias.rp.nai <- result[,16]-result[,35]
               bias.rp.asy <- result[,21]-result[,35]
               bias.cox <- result[,26]-result[,35]
               
               mb.cox <- mean(bias.cox, na.rm=TRUE)
               vb.cox<- var(bias.cox, na.rm=TRUE)
          }
          else if (ncol(result)==30){
               reject.ph <- sum(result[,28]<0.05 & !is.na(result[,28]))/sum(!is.na(result[,28]))
               bias.psu.sand <- result[,1]-result[,30]
               bias.psu.ajk <- result[,6]-result[,30]
               bias.km <- result[,11]-result[,30]
               bias.rp.nai <- result[,16]-result[,30]
               bias.rp.asy <- result[,21]-result[,30]
          }
          
          mb.psu.sand <- mean(bias.psu.sand, na.rm=TRUE)
          vb.psu.sand <- var(bias.psu.sand, na.rm=TRUE)
          
          mb.psu.ajk <- mean(bias.psu.ajk, na.rm=TRUE)
          vb.psu.ajk <- var(bias.psu.ajk, na.rm=TRUE)
          
          mb.km <- mean(bias.km, na.rm=TRUE)
          vb.km <- var(bias.km, na.rm=TRUE)
          
          mb.rp.nai <- mean(bias.rp.nai, na.rm=TRUE)
          vb.rp.nai <- var(bias.rp.nai, na.rm=TRUE)
          
          mb.rp.asy <- mean(bias.rp.asy, na.rm=TRUE)
          vb.rp.asy <- var(bias.rp.asy, na.rm=TRUE)
          
          
          if (ncol(result)==30){
               fin <- c("psu.sand.diff"=res.mean[1], "psu.sand.se"=res.mean[2],"psu.sand.lb"=res.mean[3],
                        "psu.sand.ub"=res.mean[4], "psu.sand.rej"=reject.psu.sand,"psu.sand.bias"=mb.psu.sand,
                        "psu.sand.vbias"=vb.psu.sand,
                        
                        "psu.ajk.diff"=res.mean[6],"psu.ajk.se"=res.mean[7],"psu.ajk.lb"=res.mean[8],
                        "psu.ajk.ub"=res.mean[9],"psu.ajk.rej"=reject.psu.ajk,"psu.ajk.bias"=mb.psu.ajk,
                        "psu.ajk.vbias"=vb.psu.ajk,
                        
                        "km.diff"=res.mean[11],"km.se"=res.mean[12],"km.lb"=res.mean[13],
                        "km.ub"=res.mean[14],"km.rej"=reject.km,"km.bias"=mb.km,"km.vbias"=vb.km,
                        
                        "rp.nai.diff"=res.mean[16],"rp.nai.se"=res.mean[17],"rp.nai.lb"=res.mean[18],
                        "rp.nai.ub"=res.mean[19],"rp.nai.rej"=reject.rp.nai,"rp.nai.bias"=mb.rp.nai,
                        "rp.nai.vbias"=vb.rp.nai,
                        
                        "rp.asy.diff"=res.mean[21],"rp.asy.se"=res.mean[22],"rp.asy.lb"=res.mean[23],
                        "rp.asy.ub"=res.mean[24],"rp.asy.rej"=reject.rp.asy,"rp.asy.bias"=mb.rp.asy,
                        "rp.asy.vbias"=vb.rp.asy,
                        
                        "duration"=res.mean[26],"enrolled"=res.mean[27], "ph.reject"=reject.ph, 
                        "tmax"=res.mean[29], "trmst"=res.mean[30])
          }
          
          else if (ncol(result)==35){
               fin <- c("psu.sand.diff"=res.mean[1], "psu.sand.se"=res.mean[2],"psu.sand.lb"=res.mean[3],
                        "psu.sand.ub"=res.mean[4], "psu.sand.rej"=reject.psu.sand,"psu.sand.bias"=mb.psu.sand,
                        "psu.sand.vbias"=vb.psu.sand,
                        
                        "psu.ajk.diff"=res.mean[6],"psu.ajk.se"=res.mean[7],"psu.ajk.lb"=res.mean[8],
                        "psu.ajk.ub"=res.mean[9],"psu.ajk.rej"=reject.psu.ajk,"psu.ajk.bias"=mb.psu.ajk,
                        "psu.ajk.vbias"=vb.psu.ajk,
                        
                        "tian.diff"=res.mean[11],"tian.se"=res.mean[12],"tian.lb"=res.mean[13],
                        "tian.ub"=res.mean[14],"tian.rej"=reject.km,"tian.bias"=mb.km,"tian.vbias"=vb.km,
                        
                        "rp.nai.diff"=res.mean[16],"rp.nai.se"=res.mean[17],"rp.nai.lb"=res.mean[18],
                        "rp.nai.ub"=res.mean[19],"rp.nai.rej"=reject.rp.nai,"rp.nai.bias"=mb.rp.nai,
                        "rp.nai.vbias"=vb.rp.nai,
                        
                        "rp.asy.diff"=res.mean[21],"rp.asy.se"=res.mean[22],"rp.asy.lb"=res.mean[23],
                        "rp.asy.ub"=res.mean[24],"rp.asy.rej"=reject.rp.asy,"rp.asy.bias"=mb.rp.asy,
                        "rp.asy.vbias"=vb.rp.asy,
                        
                        "cox.asy.diff"=res.mean[26],"cox.asy.se"=res.mean[27],"cox.asy.lb"=res.mean[28],
                        "cox.asy.ub"=res.mean[29],"cox.asy.rej"=reject.cox, "cox.bias"=mb.cox, "cox.vbias"=vb.cox,
                        
                        "duration"=res.mean[31],"enrolled"=res.mean[32], "ph.reject"=reject.ph,
                        "tmax"=res.mean[34], "trmst"=res.mean[35])
          }
          
          return(fin)
}



#### for sim1 nocov, 550; top half of table for scenario 1
sim1nocov <- readRDS("Results/sim1/PH.ss550trteff3nocov.nocoef.LTFU0.02.RDS")
result1top <- sum.func(sim1nocov)
result1top
#### the RP method failed to converge in 10 trials (1%)
summary(sim1nocov)
# type 1 error row
sim1type <- readRDS("Results/sim1type1/PH.ss550trteff0nocov.nocoef.LTFU0.02.RDS")
result1t1 <- sum.func(sim1type)
result1t1


### for sim1 including age and female, also 550; bottom half of scenario 1
sim1cov <- readRDS("Results/sim1/PH.ss550trteff3nocov.age.fem.LTFU0.02.RDS")
result1bot <- sum.func(sim1cov)
result1bot
#### the RP method failed to converge in 10 trials (1%)
summary(sim1cov)
#type 1 error row
sim1type1cov <- readRDS("Results/sim1type1/PH.ss550trteff0nocov.age.fem.LTFU0.02.RDS")
result1b1 <- sum.func(sim1type1cov) 
result1b1


##### Table 1 in the manuscript
tab1 <- matrix(nrow=10,ncol=7)
tab1[1,1] <- result1top["km.bias"]
tab1[2,1] <- result1top["km.vbias"]
tab1[3,1] <- result1top["km.se"]
tab1[4,1] <- result1top["km.rej"]
tab1[5,1] <- result1t1["km.rej"]

tab1[1,4] <- result1top["psu.ajk.bias"]
tab1[2,4] <- result1top["psu.ajk.vbias"]
tab1[3,4] <- result1top["psu.ajk.se"]
tab1[4,4] <- result1top["psu.ajk.rej"]
tab1[5,4] <- result1t1["psu.ajk.rej"]

tab1[1,5] <- result1top["psu.sand.bias"]
tab1[2,5] <- result1top["psu.sand.vbias"]
tab1[3,5] <- result1top["psu.sand.se"]
tab1[4,5] <- result1top["psu.sand.rej"]
tab1[5,5] <- result1t1["psu.sand.rej"]

tab1[1,6] <- result1top["rp.nai.bias"]
tab1[2,6] <- result1top["rp.nai.vbias"]
tab1[3,6] <- result1top["rp.nai.se"]
tab1[4,6] <- result1top["rp.nai.rej"]
tab1[5,6] <- result1t1["rp.nai.rej"]

tab1[1,7] <- result1top["rp.asy.bias"]
tab1[2,7] <- result1top["rp.asy.vbias"]
tab1[3,7] <- result1top["rp.asy.se"]
tab1[4,7] <- result1top["rp.asy.rej"]
tab1[5,7] <- result1t1["rp.asy.rej"]

tab1[6,2] <- result1bot["tian.bias"]
tab1[7,2] <- result1bot["tian.vbias"]
tab1[8,2] <- result1bot["tian.se"]
tab1[9,2] <- result1bot["tian.rej"]
tab1[10,2] <- result1b1["tian.rej"]

tab1[6,3] <- result1bot["cox.bias"]
tab1[7,3] <- result1bot["cox.vbias"]
tab1[8,3] <- result1bot["cox.asy.se"]
tab1[9,3] <- result1bot["cox.asy.rej"]
tab1[10,3] <- result1b1["cox.asy.rej"]

tab1[6,4] <- result1bot["psu.ajk.bias"]
tab1[7,4] <- result1bot["psu.ajk.vbias"]
tab1[8,4] <- result1bot["psu.ajk.se"]
tab1[9,4] <- result1bot["psu.ajk.rej"]
tab1[10,4] <- result1b1["psu.ajk.rej"]

tab1[6,5] <- result1bot["psu.sand.bias"]
tab1[7,5] <- result1bot["psu.sand.vbias"]
tab1[8,5] <- result1bot["psu.sand.se"]
tab1[9,5] <- result1bot["psu.sand.rej"]
tab1[10,5] <- result1b1["psu.sand.rej"]

tab1[6,6] <- result1bot["rp.nai.bias"]
tab1[7,6] <- result1bot["rp.nai.vbias"]
tab1[8,6] <- result1bot["rp.nai.se"]
tab1[9,6] <- result1bot["rp.nai.rej"]
tab1[10,6] <- result1b1["rp.nai.rej"]

tab1[6,7] <- result1bot["rp.asy.bias"]
tab1[7,7] <- result1bot["rp.asy.vbias"]
tab1[8,7] <- result1bot["rp.asy.se"]
tab1[9,7] <- result1bot["rp.asy.rej"]
tab1[10,7] <- result1b1["rp.asy.rej"]

tab1 <- round(tab1, digits=3)

tab1 <- data.frame(tab1)

colnames(tab1) <- c("KM","Tian","Chen","Psu.ajk","Psu.sand","FPM.delta","FPM.Mest")
rownames(tab1) <- c("Bias","var(Bias)","Mean(SE)","Power","Type 1 error","Bias2","var(Bias)2","Mean(SE)2","Power2","Type 1 error2")
print(tab1)

#### for sim2 nocov, ss 150; top half of table for scenario 2
sim2nocov <- readRDS("Results/sim2/PH.ss150trteff3nocov.age.neg2.fem3.LTFU0.02.RDS")
result2top <- sum.func(sim2nocov)
result2top
### All methods converged on all trials
summary(sim2nocov)
# type 1 error row
sim2type1 <- readRDS("Results/sim2type1/PH.ss150trteff0nocov.age.neg2.fem3.LTFU0.02.RDS")
result2t1 <- sum.func(sim2type1)
result2t1


#### for sim2 nocov, ss 150; bottom half of table for scenario 2
sim2cov <- readRDS("Results/sim2/PH.ss150trteff3covsall.age.neg2.fem3.LTFU0.02.RDS")
result2bot <- sum.func(sim2cov)
result2bot
### RP method failed to converge on 2 trials (.2%)
summary(sim2cov)
# type 1 error row
sim2type1cov <- readRDS("Results/sim2type1/PH.ss150trteff0covsall.age.neg2.fem3.LTFU0.02.RDS")
result2b1 <- sum.func(sim2type1cov)
result2b1


##### Table 2 in the manuscript
tab2 <- matrix(nrow=10,ncol=7)
tab2[1,1] <- result2top["km.bias"]
tab2[2,1] <- result2top["km.vbias"]
tab2[3,1] <- result2top["km.se"]
tab2[4,1] <- result2top["km.rej"]
tab2[5,1] <- result2t1["km.rej"]

tab2[1,4] <- result2top["psu.ajk.bias"]
tab2[2,4] <- result2top["psu.ajk.vbias"]
tab2[3,4] <- result2top["psu.ajk.se"]
tab2[4,4] <- result2top["psu.ajk.rej"]
tab2[5,4] <- result2t1["psu.ajk.rej"]

tab2[1,5] <- result2top["psu.sand.bias"]
tab2[2,5] <- result2top["psu.sand.vbias"]
tab2[3,5] <- result2top["psu.sand.se"]
tab2[4,5] <- result2top["psu.sand.rej"]
tab2[5,5] <- result2t1["psu.sand.rej"]

tab2[1,6] <- result2top["rp.nai.bias"]
tab2[2,6] <- result2top["rp.nai.vbias"]
tab2[3,6] <- result2top["rp.nai.se"]
tab2[4,6] <- result2top["rp.nai.rej"]
tab2[5,6] <- result2t1["rp.nai.rej"]

tab2[1,7] <- result2top["rp.asy.bias"]
tab2[2,7] <- result2top["rp.asy.vbias"]
tab2[3,7] <- result2top["rp.asy.se"]
tab2[4,7] <- result2top["rp.asy.rej"]
tab2[5,7] <- result2t1["rp.asy.rej"]

tab2[6,2] <- result2bot["tian.bias"]
tab2[7,2] <- result2bot["tian.vbias"]
tab2[8,2] <- result2bot["tian.se"]
tab2[9,2] <- result2bot["tian.rej"]
tab2[10,2] <- result2b1["tian.rej"]

tab2[6,3] <- result2bot["cox.bias"]
tab2[7,3] <- result2bot["cox.vbias"]
tab2[8,3] <- result2bot["cox.asy.se"]
tab2[9,3] <- result2bot["cox.asy.rej"]
tab2[10,3] <- result2b1["cox.asy.rej"]

tab2[6,4] <- result2bot["psu.ajk.bias"]
tab2[7,4] <- result2bot["psu.ajk.vbias"]
tab2[8,4] <- result2bot["psu.ajk.se"]
tab2[9,4] <- result2bot["psu.ajk.rej"]
tab2[10,4] <- result2b1["psu.ajk.rej"]

tab2[6,5] <- result2bot["psu.sand.bias"]
tab2[7,5] <- result2bot["psu.sand.vbias"]
tab2[8,5] <- result2bot["psu.sand.se"]
tab2[9,5] <- result2bot["psu.sand.rej"]
tab2[10,5] <- result2b1["psu.sand.rej"]

tab2[6,6] <- result2bot["rp.nai.bias"]
tab2[7,6] <- result2bot["rp.nai.vbias"]
tab2[8,6] <- result2bot["rp.nai.se"]
tab2[9,6] <- result2bot["rp.nai.rej"]
tab2[10,6] <- result2b1["rp.nai.rej"]

tab2[6,7] <- result2bot["rp.asy.bias"]
tab2[7,7] <- result2bot["rp.asy.vbias"]
tab2[8,7] <- result2bot["rp.asy.se"]
tab2[9,7] <- result2bot["rp.asy.rej"]
tab2[10,7] <- result2b1["rp.asy.rej"]

tab2 <- round(tab2, digits=3)

tab2 <- data.frame(tab2)

colnames(tab2) <- c("KM","Tian","Chen","Psu.ajk","Psu.sand","FPM.delta","FPM.Mest")
rownames(tab2) <- c("Bias","var(Bias)","Mean(SE)","Power","Type 1 error","Bias2","var(Bias)2","Mean(SE)2","Power2","Type 1 error2")
print(tab2)

