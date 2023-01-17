#### Simulation 2 comparing RMST estimation methods


##### RMSTdiff comparison simulation 1
##### this simulation explores power at different sample sizes and same treatment effect and data generated without covariates

source("Programs/sourceRMSTdiff.R")
#### checking type 1 error

reps <- 1000
ss <- 550
LTFU <- .02
b0 <- 9
bstar <- 3


accru <- 10

Sys.time()



for (k in 1:length(ss)){

     result <- matrix(nrow=reps, ncol=30)
     set.seed(12)
     for(i in 1:reps){
          
          data <- simdat.cov(n=(5*ss[k]),               #total number of observations
                             nevents=ss[k],
                             accru=accru,           # accrual rate per month
                             LTFU=LTFU,            # loss to follow up proportion
                             pos=0.5,             # proportion of biomarker positive subjects
                             trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                             female=0.5,          # probability of being female
                             age.min=55,         # minimum age of subjects
                             age.max=85,         # maximum age of subjects
                             dis1=.1,            # probability of having disease severity 1
                             dis2=.2,            # probability of having disease severity 2
                             dis3=.3,            # probability of having disease severity 3
                             dis4=.2,            # probability of having disease severity 4
                             b0=b0,              # the median survival for a male of minimum age with disease status 1
                             bstar=bstar,            # increase in median survival for being treated with B
                             bmark=0,             # coefficient for being biomarker positive
                             b1=0,              # the coefficient for age in the model to generate survival times
                             b2=0,                 # coefficient for being female
                             b3=0,              # coefficient for having disease status 2
                             b4=0,              # coefficient for having disease status 3
                             b5=0,              # coefficient for having disease status 4
                             b6=0,               # coefficient for having disease status 5
                             b7=0,               # coefficent for interaction between trt B and biomarker positive
                             PH=TRUE,            # Proportional hazards true or false; if true, both arms have exp distr
                             shapeB=3,            # shape of weibull for trt B if PH=FALSE
                             clin.trial=TRUE
          )
          
          covs <- NULL
          
          ## calculating tmax
          data <- data.frame(data)
          trt1 <- data[data$trtn==1,]
          trt0 <- data[data$trtn==0,]
          max1 <- max(trt1$surv.time)
          max0 <- max(trt0$surv.time)
          tmax <- min(max0,max1)
          
          pseudo.sand <- NULL
          pseudo.sand <- RMSTdiff(data=data,
                                  time.var="surv.time",
                                  event.var="event",
                                  trtn.var="trtn",
                                  covs=covs,
                                  tmax=tmax,
                                  method="pseudo",
                                  var.type1="sandwich",
                                  var.type2="asymp",
                                  alpha=.05)
          
          pseudo.ajk <- NULL
          pseudo.ajk <- RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="pseudo",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05)
          
          km <- NULL
          km <- RMSTdiff(data=data,
                         time.var="surv.time",
                         event.var="event",
                         trtn.var="trtn",
                         covs=covs,
                         tmax=tmax,
                         method="KM",
                         var.type1="ajk",
                         var.type2="asymp",
                         alpha=.05)
          
          rp31.naive <- NULL
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
          
          rp31.asymp <- NULL
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
          
          ## calculating the p-value from the coxph fit
          formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
          temp <- NULL
          data <- data.frame(data)
          temp <- coxph(formula,data=data)
          p <- summary(temp)$coefficients["trtn","Pr(>|z|)"]
          
          ## calculating the true RMST at the follow-up time
          lam0 <- log(2)/b0
          exp.rmst0 <- function(x){
               exp(-lam0*x)
          }
          lam1 <- log(2)/(b0+bstar)
          exp.rmst1 <- function(x){
               exp(-lam1*x)
          }
          rmst0 <- integrate(exp.rmst0, lower=0, upper=tmax)$value
          rmst1 <- integrate(exp.rmst1, lower=0, upper=tmax)$value
          trmst <- rmst1-rmst0
          
          
          result[i,] <- c(pseudo.sand, pseudo.ajk, km, rp31.naive, rp31.asymp, data[1,"dur"],data[1,"enr"],p, tmax, trmst)
          
          
     }
     
     filetmp <- paste("PH.ss",ss[k],"trteff",bstar,"nocov.nocoef.LTFU",LTFU,".RDS",sep="")
     
     saveRDS(result, file=paste("Results","sim1", filetmp, sep="/"))
     
}




for (k in 1:length(ss)){
     
     
     result <- matrix(nrow=reps, ncol=35)
     set.seed(12)
     for(i in 1:reps){
          
          data <- simdat.cov(n=(5*ss[k]),               #total number of observations
                             nevents=ss[k],
                             accru=accru,           # accrual rate per month
                             LTFU=LTFU,            # loss to follow up proportion
                             pos=0.5,             # proportion of biomarker positive subjects
                             trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                             female=0.5,          # probability of being female
                             age.min=55,         # minimum age of subjects
                             age.max=85,         # maximum age of subjects
                             dis1=.1,            # probability of having disease severity 1
                             dis2=.2,            # probability of having disease severity 2
                             dis3=.3,            # probability of having disease severity 3
                             dis4=.2,            # probability of having disease severity 4
                             b0=9,              # the median survival for a male of minimum age with disease status 1
                             bstar=bstar,            # increase in median survival for being treated with B
                             bmark=0,             # coefficient for being biomarker positive
                             b1=0,              # the coefficient for age in the model to generate survival times
                             b2=0,                 # coefficient for being female
                             b3=0,              # coefficient for having disease status 2
                             b4=0,              # coefficient for having disease status 3
                             b5=0,              # coefficient for having disease status 4
                             b6=0,               # coefficient for having disease status 5
                             b7=0,               # coefficent for interaction between trt B and biomarker positive
                             PH=TRUE,            # Proportional hazards true or false; if true, both arms have exp distr
                             shapeB=3,            # shape of weibull for trt B if PH=FALSE
                             clin.trial=TRUE
          )
          
          covs <- c("age.cent","female")
          
          ## calculating tmax
          data <- data.frame(data)
          trt1 <- data[data$trtn==1,]
          trt0 <- data[data$trtn==0,]
          max1 <- max(trt1$surv.time)
          max0 <- max(trt0$surv.time)
          tmax <- min(max0,max1)
          
          pseudo.sand <- RMSTdiff(data=data,
                                  time.var="surv.time",
                                  event.var="event",
                                  trtn.var="trtn",
                                  covs=covs,
                                  tmax=tmax,
                                  method="pseudo",
                                  var.type1="sandwich",
                                  var.type2="asymp",
                                  alpha=.05)
          
          pseudo.ajk <- RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="pseudo",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05)
          
          Tian <- RMSTdiff(data=data,
                           time.var="surv.time",
                           event.var="event",
                           trtn.var="trtn",
                           covs=covs,
                           tmax=tmax,
                           method="Tian",
                           var.type1="ajk",
                           var.type2="asymp",
                           alpha=.05)
          
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
          
          cox <- RMSTdiff(data=data,
                          time.var="surv.time",
                          event.var="event",
                          trtn.var="trtn",
                          covs=covs,
                          tmax=tmax,
                          method="Chen",
                          var.type1="ajk",
                          var.type2="asymp",
                          alpha=.05)
          
          formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
          temp <- NULL
          data <- data.frame(data)
          temp <- coxph(formula,data=data)
          p <- summary(temp)$coefficients["trtn","Pr(>|z|)"]
          
          ## calculating the true RMST at the follow-up time
          lam0 <- log(2)/b0
          exp.rmst0 <- function(x){
               exp(-lam0*x)
          }
          lam1 <- log(2)/(b0+bstar)
          exp.rmst1 <- function(x){
               exp(-lam1*x)
          }
          rmst0 <- integrate(exp.rmst0, lower=0, upper=tmax)$value
          rmst1 <- integrate(exp.rmst1, lower=0, upper=tmax)$value
          trmst <- rmst1-rmst0
          
          result[i,] <- c(pseudo.sand, pseudo.ajk, Tian, rp31.naive, rp31.asymp, cox, data[1,"dur"],data[1,"enr"],p, tmax, trmst)
          
          
     }
     
     filetmp <- paste("PH.ss",ss[k],"trteff",bstar,"nocov.age.fem.LTFU",LTFU,".RDS",sep="")
     
     saveRDS(result, file=paste("Results","sim1", filetmp, sep="/"))
     
}





for (k in 1:length(ss)){
     
     
     result <- matrix(nrow=reps, ncol=35)
     set.seed(12)
     for(i in 1:reps){
          
          data <- simdat.cov(n=(5*ss[k]),               #total number of observations
                             nevents=ss[k],
                             accru=accru,           # accrual rate per month
                             LTFU=LTFU,            # loss to follow up proportion
                             pos=0.5,             # proportion of biomarker positive subjects
                             trtB=0.5,            # proportion of subjects randomly assigned to treatment B
                             female=0.5,          # probability of being female
                             age.min=55,         # minimum age of subjects
                             age.max=85,         # maximum age of subjects
                             dis1=.1,            # probability of having disease severity 1
                             dis2=.2,            # probability of having disease severity 2
                             dis3=.3,            # probability of having disease severity 3
                             dis4=.2,            # probability of having disease severity 4
                             b0=9,              # the median survival for a male of minimum age with disease status 1
                             bstar=bstar,            # increase in median survival for being treated with B
                             bmark=0,             # coefficient for being biomarker positive
                             b1=0,              # the coefficient for age in the model to generate survival times
                             b2=0,                 # coefficient for being female
                             b3=0,              # coefficient for having disease status 2
                             b4=0,              # coefficient for having disease status 3
                             b5=0,              # coefficient for having disease status 4
                             b6=0,               # coefficient for having disease status 5
                             b7=0,               # coefficent for interaction between trt B and biomarker positive
                             PH=TRUE,            # Proportional hazards true or false; if true, both arms have exp distr
                             shapeB=3,            # shape of weibull for trt B if PH=FALSE
                             clin.trial=TRUE
          )
          
          covs <- c("female")
          
          ## calculating tmax
          data <- data.frame(data)
          trt1 <- data[data$trtn==1,]
          trt0 <- data[data$trtn==0,]
          max1 <- max(trt1$surv.time)
          max0 <- max(trt0$surv.time)
          tmax <- min(max0,max1)
          
          pseudo.sand <- RMSTdiff(data=data,
                                  time.var="surv.time",
                                  event.var="event",
                                  trtn.var="trtn",
                                  covs=covs,
                                  tmax=tmax,
                                  method="pseudo",
                                  var.type1="sandwich",
                                  var.type2="asymp",
                                  alpha=.05)
          
          pseudo.ajk <- RMSTdiff(data=data,
                                 time.var="surv.time",
                                 event.var="event",
                                 trtn.var="trtn",
                                 covs=covs,
                                 tmax=tmax,
                                 method="pseudo",
                                 var.type1="ajk",
                                 var.type2="asymp",
                                 alpha=.05)
          
          Tian <- RMSTdiff(data=data,
                         time.var="surv.time",
                         event.var="event",
                         trtn.var="trtn",
                         covs=covs,
                         tmax=tmax,
                         method="Tian",
                         var.type1="ajk",
                         var.type2="asymp",
                         alpha=.05)
          
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
          
          cox <- RMSTdiff(data=data,
                          time.var="surv.time",
                          event.var="event",
                          trtn.var="trtn",
                          covs=covs,
                          tmax=tmax,
                          method="Chen",
                          var.type1="ajk",
                          var.type2="asymp",
                          alpha=.05)
          
          formula <- as.formula(paste0("Surv(surv.time,event) ~", paste(covs, collapse = "+"),"+ trtn"))
          temp <- NULL
          data <- data.frame(data)
          temp <- coxph(formula,data=data)
          p <- summary(temp)$coefficients["trtn","Pr(>|z|)"]
          
          ## calculating the true RMST at the follow-up time
          lam0 <- log(2)/b0
          exp.rmst0 <- function(x){
               exp(-lam0*x)
          }
          lam1 <- log(2)/(b0+bstar)
          exp.rmst1 <- function(x){
               exp(-lam1*x)
          }
          rmst0 <- integrate(exp.rmst0, lower=0, upper=tmax)$value
          rmst1 <- integrate(exp.rmst1, lower=0, upper=tmax)$value
          trmst <- rmst1-rmst0
          
          result[i,] <- c(pseudo.sand, pseudo.ajk, Tian, rp31.naive, rp31.asymp, cox, data[1,"dur"],data[1,"enr"],p, tmax, trmst)
          
          
     }
     
     filetmp <- paste("PH.ss",ss[k],"trteff",bstar,"nocov.fem.LTFU",LTFU,".RDS",sep="")
     
     saveRDS(result, file=paste("Results","sim1", filetmp, sep="/"))
     
}






     