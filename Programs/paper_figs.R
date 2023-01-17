####### code to produce the figure for the paper Implementing RMST in R

source("Programs/sourceRMSTdiff.R")


### the following is the example dataset , figure and output for Section 3
set.seed(12)
data <- simdat.cov(n=2000,               #total number of observations
                   nevents=400,
                   accru=10,           # accrual rate per month
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
                   b0=12,              # the median survival for a male of minimum age with disease status 1
                   bstar=4,            # increase in median survival for being treated with B
                   bmark=0,             # coefficient for being biomarker positive
                   b1=-2,              # the coefficient for age in the model to generate survival times
                   b2=3,                 # coefficient for being female
                   b3=0,              # coefficient for having disease status 2
                   b4=0,              # coefficient for having disease status 3
                   b5=0,              # coefficient for having disease status 4
                   b6=0,               # coefficient for having disease status 5
                   b7=0,               # coefficent for interaction between trt B and biomarker positive
                   PH=FALSE,            # Proportional hazards true or false; if true, both arms have exp distr
                   shapeB=3,            # shape of weibull for trt B if PH=FALSE
                   clin.trial=FALSE
)


data <- data.frame(data)

data <- data[,c("id","entry.time","trtn","female","age.cent","event","surv.time","event.time")]
sink(file="output1.txt")
head(data)
closeAllConnections()

sink(file="output2.txt")
RMSTdiff(data=data,
         time.var="surv.time",
         event.var="event",
         trtn.var="trtn",
         covs=c("female","age.cent"),
         tmax=20,
         method="pseudo",
         var.type1="ajk")
closeAllConnections()


Y <- Surv(data$surv.time, data$event)
kmfit <- survfit(Y~data$trtn)
plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,30), xlab="Time in Months",
     ylab="Survival Probability")
abline(v=20,col="red")
legend(22,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)


sink(file="output3.txt")
undebug(RMSTdiff)
RMSTdiff(data=data,
         time.var="surv.time",
         event.var="event",
         trtn.var="trtn",
         covs=c("female","age.cent"),
         tmax=20,
         method="pseudo",
         var.type1="sandwich")
closeAllConnections()


sink(file="output4.txt")
RMSTdiff(data=data,
         time.var="surv.time",
         event.var="event",
         trtn.var="trtn",
         covs=c("female","age.cent"),
         tmax=20,
         method="Cox",
         var.type1="sandwich",
         var.type2="asymp")
closeAllConnections()

### figure 1
x <- seq(0,60, by=.01)
y <- exp(-(log(2)/18)*x)
z <- exp(-((log(2)/18)^3)*x^3)

plot(x,y,type="l", col="blue", xlab = "Survival Time In Months", 
     ylab = "Survival Probability", main="Difference in Restricted Mean Survival at 24 Months")
lines(x, z,type="l", col="purple")
abline(v=24, col="red")
polygon(c(x[x<=24],rev(x[x<=24])),c(z[x<=24],rev(y[x<=24])),col="green")
legend(40,0.9,legend=c("Group 1", "Group 0"), col=c("purple","blue"),lty=1:1)



### figure for scenario 1
curve(exp(-(log(2)/12)*x), from=0,to=60, col="green",lty=2, xlab="Time in Months",
      ylab="Survival Probability", main="Scenario 1")
curve(exp(-(log(2)/9)*x), from=0,to=60, col="blue",lty=1, add=TRUE)
legend(30,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)


# ### figure for scenario 2
# 
# set.seed(1212)
# data2 <- simdat.cov(n=1000000,               #total number of observations
#                     nevents=400,
#                     accru=10,           # accrual rate per month
#                     LTFU=0,            # loss to follow up proportion
#                     pos=0.5,             # proportion of biomarker positive subjects
#                     trtB=0.5,            # proportion of subjects randomly assigned to treatment B
#                     female=0.5,          # probability of being female
#                     age.min=55,         # minimum age of subjects
#                     age.max=85,         # maximum age of subjects
#                     dis1=.2,            # probability of having disease severity 1
#                     dis2=.2,            # probability of having disease severity 2
#                     dis3=.3,            # probability of having disease severity 3
#                     dis4=.2,            # probability of having disease severity 4
#                     b0=12,              # the median survival for a male of minimum age with disease status 1
#                     bstar=9,            # increase in median survival for being treated with B
#                     bmark=0,             # coefficient for being biomarker positive
#                     b1=0,              # the coefficient for age in the model to generate survival times
#                     b2=0,                 # coefficient for being female
#                     b3=0,              # coefficient for having disease status 2
#                     b4=0,              # coefficient for having disease status 3
#                     b5=0,              # coefficient for having disease status 4
#                     b6=0,               # coefficient for having disease status 5
#                     b7=0,               # coefficent for interaction between trt B and biomarker positive
#                     PH=FALSE,            # Proportional hazards true or false; if true, both arms have exp distr
#                     shapeB=3,            # shape of weibull for trt B if PH=FALSE
#                     clin.trial=FALSE
# )
# 
# data <- data.frame(data2)
# Y <- Surv(data$surv.time, data$event)
# 
# kmfit <- survfit(Y~data$trtn)
# 
# plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,60),
#      xlab="Time in Months",ylab="Survival Probability",main="Scenario 2")
# legend(30,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)
# 

### KM curve for scenario 2
set.seed(12)
data2 <- simdat.cov(n=1000000,               #total number of observations
                    nevents=400,
                    accru=10,           # accrual rate per month
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
                    bstar=3,            # increase in median survival for being treated with B
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

data <- data.frame(data2)
Y <- Surv(data$surv.time, data$event)

kmfit <- survfit(Y~data$trtn)

plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,60),
     xlab="Time in Months",ylab="Survival Probability",main="Scenario 2")
legend(30,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)



### KM curve for scenario 4
# set.seed(1212)
# data2 <- simdat.cov(n=1000000,               #total number of observations
#                     nevents=400,
#                     accru=10,           # accrual rate per month
#                     LTFU=0,            # loss to follow up proportion
#                     pos=0.5,             # proportion of biomarker positive subjects
#                     trtB=0.5,            # proportion of subjects randomly assigned to treatment B
#                     female=0.5,          # probability of being female
#                     age.min=55,         # minimum age of subjects
#                     age.max=85,         # maximum age of subjects
#                     dis1=.2,            # probability of having disease severity 1
#                     dis2=.2,            # probability of having disease severity 2
#                     dis3=.3,            # probability of having disease severity 3
#                     dis4=.2,            # probability of having disease severity 4
#                     b0=12,              # the median survival for a male of minimum age with disease status 1
#                     bstar=7,            # increase in median survival for being treated with B
#                     bmark=0,             # coefficient for being biomarker positive
#                     b1=-2,              # the coefficient for age in the model to generate survival times
#                     b2=3,                 # coefficient for being female
#                     b3=0,              # coefficient for having disease status 2
#                     b4=0,              # coefficient for having disease status 3
#                     b5=0,              # coefficient for having disease status 4
#                     b6=0,               # coefficient for having disease status 5
#                     b7=0,               # coefficent for interaction between trt B and biomarker positive
#                     PH=FALSE,            # Proportional hazards true or false; if true, both arms have exp distr
#                     shapeB=3,            # shape of weibull for trt B if PH=FALSE
#                     clin.trial=FALSE
# )
# 
# data <- data.frame(data2)
# Y <- Surv(data$surv.time, data$event)
# 
# kmfit <- survfit(Y~data$trtn)
# 
# plot(kmfit, lty=c("solid","dashed"), col=c("blue","green"),xlim=c(0,60),
#      xlab="Time in Months",ylab="Survival Probability",main="Scenario 4")
# legend(30,0.9,legend=c("Group 1", "Group 0"), col=c("green","blue"),lty=2:1)








