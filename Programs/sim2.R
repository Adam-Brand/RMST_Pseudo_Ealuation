#### Simulation 2 comparing RMST estimation methods


##### RMSTdiff comparison simulation 2
##### this simulation explores power at different sample sizes and same treatment effect and data generated with covariates

source("Programs/sourceRMSTdiff.R")
#### checking type 1 error

reps <- 1000
ss <- 150
LTFU <- .02
bstar <- 3

accru <- 10

Sys.time()

set.seed(12)
data2 <- simdat.cov(n=10000000,               #total number of observations
                   nevents=ss,
                   accru=accru,           # accrual rate per month
                   LTFU=LTFU,            # loss to follow up proportion
                   pos=0.5,             # proportion of biomarker positive subjects
                   trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                   female=0.5,          # probability of being female
                   age.min=55,         # minimum age of subjects
                   age.max=85,         # maximum age of subjects
                   dis1=.2,            # probability of having disease severity 1
                   dis2=.2,            # probability of having disease severity 2
                   dis3=.3,            # probability of having disease severity 3
                   dis4=.2,            # probability of having disease severity 4
                   b0=9,              # the median survival for a male of minimum age with disease status 1
                   bstar=bstar,            # increase in median survival for being treated with B
                   bmark=0,             # coefficient for being biomarker positive
                   b1=-2,              # the coefficient for age in the model to generate survival times
                   b2=3,                 # coefficient for being female
                   b3=0,              # coefficient for having disease status 2
                   b4=0,              # coefficient for having disease status 3
                   b5=0,              # coefficient for having disease status 4
                   b6=0,               # coefficient for having disease status 5
                   b7=0,               # coefficent for interaction between trt B and biomarker positive
                   PH=TRUE,            # Proportional hazards true or false; if true, both arms have exp distr
                   shapeB=3,            # shape of weibull for trt B if PH=FALSE
                   clin.trial=FALSE
)

data2 <- data2[order(data2[,"entry.time"]),]

data3 <- simdat.cov(n=10000000,               #total number of observations
                    nevents=ss,
                    accru=accru,           # accrual rate per month
                    LTFU=0,            # loss to follow up proportion
                    pos=0.5,             # proportion of biomarker positive subjects
                    trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                    female=0.5,          # probability of being female
                    age.min=55,         # minimum age of subjects
                    age.max=85,         # maximum age of subjects
                    dis1=.2,            # probability of having disease severity 1
                    dis2=.2,            # probability of having disease severity 2
                    dis3=.3,            # probability of having disease severity 3
                    dis4=.2,            # probability of having disease severity 4
                    b0=9,              # the median survival for a male of minimum age with disease status 1
                    bstar=bstar,            # increase in median survival for being treated with B
                    bmark=0,             # coefficient for being biomarker positive
                    b1=-2,              # the coefficient for age in the model to generate survival times
                    b2=3,                 # coefficient for being female
                    b3=0,              # coefficient for having disease status 2
                    b4=0,              # coefficient for having disease status 3
                    b5=0,              # coefficient for having disease status 4
                    b6=0,               # coefficient for having disease status 5
                    b7=0,               # coefficent for interaction between trt B and biomarker positive
                    PH=TRUE,            # Proportional hazards true or false; if true, both arms have exp distr
                    shapeB=3,            # shape of weibull for trt B if PH=FALSE
                    clin.trial=FALSE
)

data3 <- data.frame(data3)


for (k in 1:length(ss)){

     result <- matrix(nrow=reps, ncol=30)
     set.seed(12)
     for(i in 1:reps){
          
          # generating time between accruing each patient based on the accrual rate
          times <- rexp(n=(5*ss[k]), rate=accru)
          # generating study entry times based on time between recruiting each patient
          entry <- cumsum(times)-times[1]
          
          ind <- sample.int(10000000, size=(5*ss[k]))
          data <- data2[ind,]
          data[,"entry.time"] <- entry
          event.t <- data[,"surv.time"]+data[,"entry.time"]
          data[,"event.time"] <- event.t
          
          temp <- data[data[,"event"]==1,]
          temp <- temp[order(temp[,"event.time"],decreasing=FALSE),]
          end.time <- as.numeric(temp[ss[k],"event.time"])
          
          data <- data[data[,"entry.time"]<=end.time,]
          
          for (j in 1:nrow(data)){
               if (data[j,"event.time"]>end.time){data[j,"event"] =0
               data[j,"surv.time"]=end.time - data[j,"entry.time"]}
          }
          
          covs <- NULL
          
          ## calculating tmax
          data <- data.frame(data)
          trt1 <- data[data$trtn==1,]
          trt0 <- data[data$trtn==0,]
          max1 <- max(trt1$surv.time)
          max0 <- max(trt0$surv.time)
          tmax <- min(max0,max1)
          
          pseudo.sand <- try(RMSTdiff(data=data,
                                  time.var="surv.time",
                                  event.var="event",
                                  trtn.var="trtn",
                                  covs=covs,
                                  tmax=tmax,
                                  method="pseudo",
                                  var.type1="sandwich",
                                  var.type2="asymp",
                                  alpha=.05))
          if(inherits(pseudo.sand, "try-error")){
               pseudo.sand <- rep(NA,5)
          }
          
          pseudo.ajk <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="pseudo",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05))
          if(inherits(pseudo.ajk, "try-error")){
               pseudo.ajk <- rep(NA,5)
          }
          
          km <- try(RMSTdiff(data=data,
                         time.var="surv.time",
                         event.var="event",
                         trtn.var="trtn",
                         covs=covs,
                         tmax=tmax,
                         method="KM",
                         var.type1="ajk",
                         var.type2="asymp",
                         alpha=.05))
          if(inherits(km, "try-error")){
               km <- rep(NA,5)
          }
          
          rp31.naive <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="RP31",
                                 var.type1="ajk",
                                 var.type2="delta",
                                 alpha=.05))
          if(inherits(rp31.naive, "try-error")){
               rp31.naive <- rep(NA,5)
          }
          
          rp31.asymp <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="RP31",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05))
          if(inherits(rp31.asymp, "try-error")){
               rp31.asymp <- rep(NA,5)
          }
          
          formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
          temp <- NULL
          data <- data.frame(data)
          temp <- coxph(formula,data=data)
          p <- summary(temp)$coefficients["trtn","Pr(>|z|)"]
          
          ## calculating the true RMST at the follow-up time
          fit <- rmst2(time=data3$surv.time, status=data3$event, arm=data3$trtn, tau=tmax, alpha=.05)
          trmst <- fit$unadjusted.result["RMST (arm=1)-(arm=0)","Est."]
          
          result[i,] <- c(pseudo.sand, pseudo.ajk, km, rp31.naive, rp31.asymp, end.time,nrow(data),p, tmax,trmst)
          
          
     }
     
     filetmp <- paste("PH.ss",ss[k],"trteff",bstar,"nocov.age.neg2.fem3.LTFU",LTFU,".RDS",sep="")
     
     saveRDS(result, file=paste("Results","sim2", filetmp, sep="/"))
     
}


for (k in 1:length(ss)){
     
     result <- matrix(nrow=reps, ncol=35)
     set.seed(12)
     for(i in 1:reps){
          
          # generating time between accruing each patient based on the accrual rate
          times <- rexp(n=(5*ss[k]), rate=accru)
          # generating study entry times based on time between recruiting each patient
          entry <- cumsum(times)-times[1]
          
          ind <- sample.int(1000000, size=(5*ss[k]))
          data <- data2[ind,]
          data[,"entry.time"] <- entry
          event.t <- data[,"surv.time"]+data[,"entry.time"]
          data[,"event.time"] <- event.t
          
          temp <- data[data[,"event"]==1,]
          temp <- temp[order(temp[,"event.time"],decreasing=FALSE),]
          end.time <- as.numeric(temp[ss[k],"event.time"])
          
          data <- data[data[,"entry.time"]<=end.time,]
          
          for (j in 1:nrow(data)){
               if (data[j,"event.time"]>end.time){data[j,"event"] =0
               data[j,"surv.time"]=end.time - data[j,"entry.time"]}
          }
          
          covs <- c("age.cent","female")
          
          ## calculating tmax
          data <- data.frame(data)
          trt1 <- data[data$trtn==1,]
          trt0 <- data[data$trtn==0,]
          max1 <- max(trt1$surv.time)
          max0 <- max(trt0$surv.time)
          tmax <- min(max0,max1)
          
          pseudo.sand <- try(RMSTdiff(data=data,
                                  time.var="surv.time",
                                  event.var="event",
                                  trtn.var="trtn",
                                  covs=covs,
                                  tmax=tmax,
                                  method="pseudo",
                                  var.type1="sandwich",
                                  var.type2="asymp",
                                  alpha=.05))
          if(inherits(pseudo.sand, "try-error")){
               pseudo.sand <- rep(NA,5)
          }
          
          pseudo.ajk <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="pseudo",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05))
          if(inherits(pseudo.ajk, "try-error")){
               pseudo.ajk <- rep(NA,5)
          }
          
          Tian <- try(RMSTdiff(data=data,
                         time.var="surv.time",
                         event.var="event",
                         trtn.var="trtn",
                         covs=covs,
                         tmax=tmax,
                         method="Tian",
                         var.type1="ajk",
                         var.type2="asymp",
                         alpha=.05))
          if(inherits(km, "try-error")){
               km <- rep(NA,5)
          }
          
          rp31.naive <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="RP31",
                                 var.type1="ajk",
                                 var.type2="delta",
                                 alpha=.05))
          if(inherits(rp31.naive, "try-error")){
               rp31.naive <- rep(NA,5)
          }
          
          rp31.asymp <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="RP31",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05))
          if(inherits(rp31.asymp, "try-error")){
               rp31.asymp <- rep(NA,5)
          }
          
          cox <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="Chen",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05))
          if(inherits(cox, "try-error")){
               cox <- rep(NA,5)
          }
          
          formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
          temp <- NULL
          data <- data.frame(data)
          temp <- coxph(formula,data=data)
          p <- summary(temp)$coefficients["trtn","Pr(>|z|)"]
          
          ## calculating the true RMST at the follow-up time
          fit <- rmst2(time=data3$surv.time, status=data3$event, arm=data3$trtn, tau=tmax, alpha=.05)
          trmst <- fit$unadjusted.result["RMST (arm=1)-(arm=0)","Est."]
          
          result[i,] <- c(pseudo.sand, pseudo.ajk, Tian, rp31.naive, rp31.asymp, cox, end.time,nrow(data),p, tmax,trmst)
          
          
     }
     
     filetmp <- paste("PH.ss",ss[k],"trteff",bstar,"covsall.age.neg2.fem3.LTFU",LTFU,".RDS",sep="")
     
     saveRDS(result, file=paste("Results","sim2", filetmp, sep="/"))
     
}



for (k in 1:length(ss)){
     
     result <- matrix(nrow=reps, ncol=35)
     set.seed(12)
     for(i in 1:reps){
          
          # generating time between accruing each patient based on the accrual rate
          times <- rexp(n=(5*ss[k]), rate=accru)
          # generating study entry times based on time between recruiting each patient
          entry <- cumsum(times)-times[1]
          
          ind <- sample.int(1000000, size=(5*ss[k]))
          data <- data2[ind,]
          data[,"entry.time"] <- entry
          event.t <- data[,"surv.time"]+data[,"entry.time"]
          data[,"event.time"] <- event.t
          
          temp <- data[data[,"event"]==1,]
          temp <- temp[order(temp[,"event.time"],decreasing=FALSE),]
          end.time <- as.numeric(temp[ss[k],"event.time"])
          
          data <- data[data[,"entry.time"]<=end.time,]
          
          for (j in 1:nrow(data)){
               if (data[j,"event.time"]>end.time){data[j,"event"] =0
               data[j,"surv.time"]=end.time - data[j,"entry.time"]}
          }
          
          covs <- c("female")
          
          ## calculating tmax
          data <- data.frame(data)
          trt1 <- data[data$trtn==1,]
          trt0 <- data[data$trtn==0,]
          max1 <- max(trt1$surv.time)
          max0 <- max(trt0$surv.time)
          tmax <- min(max0,max1)
          
          pseudo.sand <- try(RMSTdiff(data=data,
                                  time.var="surv.time",
                                  event.var="event",
                                  trtn.var="trtn",
                                  covs=covs,
                                  tmax=tmax,
                                  method="pseudo",
                                  var.type1="sandwich",
                                  var.type2="asymp",
                                  alpha=.05))
          if(inherits(pseudo.sand, "try-error")){
               pseudo.sand <- rep(NA,5)
          }
          
          pseudo.ajk <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="pseudo",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05))
          if(inherits(pseudo.ajk, "try-error")){
               pseudo.ajk <- rep(NA,5)
          }
          
          Tian <- try(RMSTdiff(data=data,
                         time.var="surv.time",
                         event.var="event",
                         trtn.var="trtn",
                         covs=covs,
                         tmax=tmax,
                         method="Tian",
                         var.type1="ajk",
                         var.type2="asymp",
                         alpha=.05))
          if(inherits(km, "try-error")){
               km <- rep(NA,5)
          }
          
          rp31.naive <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="RP31",
                                 var.type1="ajk",
                                 var.type2="delta",
                                 alpha=.05))
          if(inherits(rp31.naive, "try-error")){
               rp31.naive <- rep(NA,5)
          }
          
          rp31.asymp <- try(RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="RP31",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05))
          if(inherits(rp31.asymp, "try-error")){
               rp31.asymp <- rep(NA,5)
          }
      
          cox <- try(RMSTdiff(data=data,
                          time.var="surv.time",
                          event.var="event",
                          trtn.var="trtn",
                          covs=covs,
                          tmax=tmax,
                          method="Chen",
                          var.type1="ajk",
                          var.type2="asymp",
                          alpha=.05))
          if(inherits(cox, "try-error")){
               cox <- rep(NA,5)
          }
          
          formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
          temp <- NULL
          data <- data.frame(data)
          temp <- coxph(formula,data=data)
          p <- summary(temp)$coefficients["trtn","Pr(>|z|)"]
          
          ## calculating the true RMST at the follow-up time
          fit <- rmst2(time=data3$surv.time, status=data3$event, arm=data3$trtn, tau=tmax, alpha=.05)
          trmst <- fit$unadjusted.result["RMST (arm=1)-(arm=0)","Est."]
          
          result[i,] <- c(pseudo.sand, pseudo.ajk, Tian, rp31.naive, rp31.asymp, cox, end.time,nrow(data),p, tmax,trmst)
          
          
     }
     
     filetmp <- paste("PH.ss",ss[k],"trteff",bstar,"covs.fem.age.neg2.fem3.LTFU",LTFU,".RDS",sep="")
     
     saveRDS(result, file=paste("Results","sim2", filetmp, sep="/"))
     
}


