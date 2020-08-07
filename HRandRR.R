### HRandRR.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Jun 29 2020 (14:15) 
## Version: 
## Last-Updated: Aug  7 2020 (13:22) 
##           By: Paul Blanche
##     Update #: 230
#----------------------------------------------------------------------
## 
### Commentary:
##
##  This R code produces all plots of the manuscript entitled "A hazard ratio
##  above one does not necessarily mean higher risk, when using a time-dependent Cox model".
##  It also produces data corresponding to the example and the estimates of 
##  risks and risk and hazard ratios. Theoretical (i.e. true) generated risks and hazard ratios
##  are also computed using closed form formulas for comparison and check.
##  The code loads 2 packages and sources 2 scripts of helper functions.
##
#----------------------------------------------------------------------
## 
### Code:

rm(list=ls())

setwd("~/research/Cox-time-dep-and-multistate/Rcode/")


# {{{ parameters
s <- 2      # landmark time
t <- 5      # prediction horizon
n <- 50000  # generated sample size
SavePlot <- FALSE
# }}}

# {{{ load helper functions
source("Generate-data.R")
source("TheoreticalProbabilities.R")
# }}}

# {{{ packages
library(prodlim)
library(survival)
# }}}

# {{{ generate example data
d <- GenerateData(n=n)
head(d)
# }}}

# {{{ Summaries of the generated data
event <- rep(1,n)
fitKM <- prodlim(Hist(time.death,event)~1,data=d)
par(mfrow=c(2,1))
plot(fitKM,xlim=quantile(d$time.death,c(0,0.9)))
abline(v=s,lty=2,col="red")
abline(v=s+t,lty=2,col="red")
barplot(table(d$trajectory))
# }}}

# {{{ Empirical conditional risks
# (i.e. simple empirical proportion observed from the genareted data)
ds <- d[d$time.death>s,]
ds$X1s <- 0
ds$X1s[which(ds$time.com1<=s)] <- 1
ds$X2s <- 0
ds$X2s[which(ds$time.com2<=s)] <- 1
## head(ds)
## summary(ds)
#
S0 <- mean(ds$time.death[ds$X1s==0 & ds$X2s==0] >s+t)
S1 <- mean(ds$time.death[ds$X1s==1 & ds$X2s==0] >s+t)
S2 <- mean(ds$time.death[ds$X1s==0 & ds$X2s==1] >s+t)
Sboth <- mean(ds$time.death[ds$X1s==1 & ds$X2s==1] >s+t)
#
cat("\n Empirical t-year risk (t=",t,"), computed at s=",s," (n=",nrow(ds),"still alive at s):\n",
    "\n When no Comorb at s:",round(100-S0*100,1),"% (=",sum(ds$time.death[ds$X1s==0 & ds$X2s==0]<=s+t),"/",sum(ds$X1s==0 & ds$X2s==0),")\n",
    "\n When comorb 1 at s only:",round(100-S1*100,1),"% (=",sum(ds$time.death[ds$X1s==1 & ds$X2s==0]<=s+t),"/",sum(ds$X1s==1 & ds$X2s==0),")\n",
    "\n When comorb 2 at s only:",round(100-S2*100,1),"% (=",sum(ds$time.death[ds$X1s==0 & ds$X2s==1]<=s+t),"/",sum(ds$X1s==0 & ds$X2s==1),")\n",
    "\n When both Comorb at s:",round(100-Sboth*100,1),"% (=",sum(ds$time.death[ds$X1s==1 & ds$X2s==1]<=s+t),"/",sum(ds$X1s==1 & ds$X2s==1),")\n",
    "\n That is, RR_1=",round((1-S1)/(1-S0),2),"and RR_2=",round((1-S2)/(1-S0),2),". \n"
    )
# }}}

# {{{ Hazard ratio estimates from a time-dependent Cox model
# prepare the data to fit the model
dcox <- tmerge(data1=d, data2=d, id=id, tstop=time.death)
dcox <- tmerge(dcox, d, id=id, X1 = tdc(time.com1))
dcox <- tmerge(dcox, d, id=id, X2 = tdc(time.com2))
dcox <- tmerge(dcox, d, id=id, event = event(time.death))
head(dcox)
# fit time-dependent Cox model
cox <- coxph(Surv(tstart, tstop, event) ~ X1+X2, data=dcox)
summary(cox)
# hazard ratio estimates
res <- cbind(exp(coef(cox)),exp(confint(cox)))
colnames(res) <- c("est","lower","upper")
res
# }}}

# {{{ check: compare observed risks vs true generated risks

# trick to facilitate the computation
d$TimeComob2 <- d$time.com2
d$TimeComob2[is.na(d$TimeComob2)] <- Inf
d$TimeComob1 <- d$time.com1
d$TimeComob1[is.na(d$TimeComob1)] <- Inf

# Initialize matrix of estimates
thenames <- c("R0",
              "R1",
              "R2",
              "R1and2",
              "P.Still.Healthy",
              "P.Healthy.To.Com1",
              "P.Healthy.To.Com2",
              "P.Healthy.To.Com1.and.Com2",
              "P.Still.Only.Com1",
              "P.Still.Only.Com2",
              "P.Com1.to.Com1.and.2",
              "P.Com2.to.Com1.and.2",
              "RR_1",
              "RR_2"
              )
RES <- matrix(NA,ncol=2,nrow=length(thenames))
rownames(RES) <- thenames
colnames(RES) <- c("Proportion","Truth (Formula)")
RES


RES["P.Still.Healthy",2] <- PStillHealthy(t)
RES["P.Still.Healthy",1] <- mean(d$time.death >t & d$TimeComob1>t & d$TimeComob2>t)

RES["P.Healthy.To.Com1",2] <-  PHealthyToCom1(t)
RES["P.Healthy.To.Com1",1] <- mean( d$time.death >t  & d$TimeComob1<t & d$TimeComob2>t )

RES["P.Healthy.To.Com2",2] <-  PHealthyToCom2(t)
RES["P.Healthy.To.Com2",1] <- mean( d$time.death >t  & d$TimeComob2<t & d$TimeComob1>t )

RES["P.Still.Only.Com1",2] <- PStillOnlyCom1(t)
RES["P.Still.Only.Com1",1] <- sum(d$time.death >s+t & d$TimeComob1<s & d$TimeComob2>(s+t))/sum(d$time.death >s & d$TimeComob1<s & d$TimeComob2>s)
    
RES["P.Com1.to.Com1.and.2",2] <- PCom1toCom1and2(t)
RES["P.Com1.to.Com1.and.2",1] <- mean(d$time.death >s+t & d$TimeComob1<s & d$TimeComob2<s+t &  d$TimeComob2>s )/mean(d$time.death >s &  d$TimeComob1<s & d$TimeComob2>s )

RES["R1",2] <- Risk1(t)
## 1-(PStillOnlyCom1(t)+PCom1toCom1and2(t)) 
RES["R1",1] <- mean(d$time.death >s & d$time.death <s+t & d$TimeComob1<s &  d$TimeComob2>s )/mean(d$time.death >s &  d$TimeComob1<s & d$TimeComob2>s )

RES["P.Still.Only.Com2",1] <- sum(d$time.death >s+t & d$TimeComob2<s & d$TimeComob1>(s+t))/sum(d$time.death >s & d$TimeComob2<s & d$TimeComob1>s)
RES["P.Still.Only.Com2",2] <- PStillOnlyCom2(t)

RES["P.Com2.to.Com1.and.2",1] <- mean(d$time.death >s+t & d$TimeComob2<s & d$TimeComob1<s+t &  d$TimeComob1>s )/mean(d$time.death >s &  d$TimeComob2<s & d$TimeComob1>s )
RES["P.Com2.to.Com1.and.2",2] <- PCom2toCom1and2(t)

## 1-(PStillOnlyCom2(t)+PCom2toCom1and2(t))  
RES["R2",1] <- mean(d$time.death >s & d$time.death <s+t & d$TimeComob2<s &  d$TimeComob1>s )/mean(d$time.death >s &  d$TimeComob2<s & d$TimeComob1>s )
RES["R2",2] <- Risk2(t)

RES["P.Healthy.To.Com1.and.Com2",2] <-  PHealthyToCom1and2(t)
RES["P.Healthy.To.Com1.and.Com2",1] <-  mean(d$time.death >(s+t)  & d$TimeComob2>s & d$TimeComob1>s & d$TimeComob2<(s+t) &  d$TimeComob1<(s+t) )/mean(d$time.death >s &  d$TimeComob2>s & d$TimeComob1>s )

RES["R1and2",2] <- Risk1and2(t)
RES["R1and2",1] <- mean(d$time.death >s & d$time.death <s+t & d$TimeComob2<s &  d$TimeComob1<s )/mean(d$time.death >s &  d$TimeComob2<s & d$TimeComob1<s )

RES["R0",] <- 1 - RES["P.Still.Healthy",] - RES["P.Healthy.To.Com1",] - RES["P.Healthy.To.Com2",] - RES["P.Healthy.To.Com1.and.Com2",]

RES["RR_1",] <- (RES["R1",]/RES["R0",]/100)
RES["RR_2",] <- (RES["R2",]/RES["R0",]/100)

# print comparisons
RES <- RES*100 
round(RES,2)
# }}}


#---------- plots ------


# {{{ Compute Risk and RR over time to plot
allttimest <- seq(0.001,t,length.out=1000)
# risks
AllRisk1 <- sapply(allttimest,Risk1)
AllRisk2 <- sapply(allttimest,Risk2)
AllRisk0 <- sapply(allttimest,function(t){1-(PStillHealthy(t) + PHealthyToCom1(t) + PHealthyToCom2(t) + PHealthyToCom1and2(t))})
AllRiskBoth <- sapply(allttimest,function(t)Risk1and2(t))
# risks ratios 
AllRR1 <- AllRisk1/AllRisk0
AllRR2 <- AllRisk2/AllRisk0
# }}}

# {{{ Curves of RR and HR over time
HR1 <- 6
HR2 <- 2
RR1 <- RES["RR_1",2]
RR2 <- RES["RR_2",2]
col1 <- "red"
col2 <- "springgreen3"
#---
if(SavePlot){pdf("RR-HR.pdf",width=8,height=7)}
par(mfrow=c(1,1),mai=c(0.9,0.9,0.4,0.4))
plot(allttimest,AllRR2,xlab="time t (years)",ylab="Ratio (RR or HR)",type="l",lwd=2,col=col2,ylim=c(0.7,6),axes=FALSE)
lines(allttimest,AllRR1,lwd=2,col=col1)
abline(v=0,col="grey",lwd=1,lty=2)
abline(h=HR2,col=col2,lwd=2,lty=2)
abline(h=HR1,col=col1,lwd=2,lty=2)
abline(h=1,col="grey",lwd=1,lty=2)
if(min(AllRR2)<1 & max(AllRR2)>1){
    xtcross <- allttimest[max(which(AllRR2>1))]
    abline(v=xtcross,col="grey",lwd=1,lty=2)
    axis(1,at=xtcross,label=round(xtcross,1))
}
legend("right",legend=c("comorbidity 1","comorbidity 2","Risk Ratio (RR)","Hazard Ratio (HR)"),
       bty="n",       
       lty=c(1,1,1,2),
       lwd=c(5,5,2,2),
       col=c(col1,col2,"black","black"))
axis(1,at=c(0,t))
axis(2,las=2,at=c(1,HR1,HR2,RR1,RR2),labels=c(1,HR1,HR2,format(RR1,nsmall=2,digits=2),format(RR2,nsmall=2,digits=2)))
if(SavePlot){dev.off()}
# }}}

# {{{ Cuves of t-year risks over time t
AddLinearApprox <- FALSE
col0 <- "blue"
col12 <- "grey"
a0 <- 0.15/(5/1.25)
#---
if(SavePlot){pdf(paste0(ifelse(AddLinearApprox,"Enriched-",""),"Risks.pdf"),width=8,height=7)}
par(mfrow=c(1,1),mai=c(0.9,0.9,0.4,0.4))
plot(allttimest,
     AllRisk1,
     xlab="time t (years)",
     ylab="t-year risk (%)",type="l",lwd=2,col=col1,ylim=c(0,1),axes=FALSE)
axis(1,at=c(0,t))
myyat <- (0:5)/5
lines(allttimest,AllRisk0,lwd=2,col=col0)
lines(allttimest,AllRisk2,lwd=2,col=col2)
lines(allttimest,AllRiskBoth,lwd=2,col=col12)
if(min(AllRR2)<1 & max(AllRR2)>1){
    xtcross <- allttimest[max(which(AllRR2>1))]
    segments(x0=xtcross,
            x1=xtcross,
            y0=-1,
            y1=AllRisk2[max(which(AllRR2>1))],
            col="grey",lwd=1,lty=2)
    segments(x0=-1,
             x1=xtcross,
             y0=AllRisk2[max(which(AllRR2>1))],
             y1=AllRisk2[max(which(AllRR2>1))],
             col="grey",lwd=1,lty=2)
    axis(1,at=xtcross,label=round(xtcross,1))
}
sapply(c(max(AllRisk2),max(AllRisk1),max(AllRisk0),max(AllRiskBoth)),
       function(h)segments(x0=-10,
                           x1=t,
                           y0=h,
                           y1=h,
                           col="grey",lwd=1,lty=2)
       )
segments(x0=t,
         x1=t,
         y0=-1,
         y1=max(c(AllRisk2,AllRisk1,AllRisk0,AllRiskBoth)),
         col="grey",lwd=1,lty=2)
myyat <- c(0,max(AllRisk2),max(AllRisk1),max(AllRisk0),max(AllRiskBoth),AllRisk2[max(which(AllRR2>1))],1)
axis(2,las=2,at=myyat,labels=round(100*myyat))
legend("topleft",
       legend=c("Both Comorbidities","Comorbidity 1","No comorbidity","Comorbidity 2"),
       col=c(col12,col1,col0,col2),
       bty="n",
       border="white",
       lwd=2,ncol=1)
if(AddLinearApprox){
    maxtsmall <- t/6
    curve(a0*x, from = 0, to = maxtsmall,lty=2,add=TRUE,lwd=3)
    curve(a0*HR1*x, from = 0, to = maxtsmall,lty=2,add=TRUE,lwd=3)
    curve(a0*HR2*x, from = 0, to = maxtsmall,lty=2,add=TRUE,lwd=3)
    curve(a0*HR2*HR1*x, from = 0, to = maxtsmall,lty=2,add=TRUE,lwd=3)
}
if(SavePlot){dev.off()}
# }}}

# {{{ Stacked plot of the proportion in each state over time

allttimest2 <- seq(0.001,3*t,length.out=100)
InC1 <- sapply(allttimest2,PHealthyToCom1)
InC2 <- sapply(allttimest2,PHealthyToCom2)
InH <- sapply(allttimest2,PStillHealthy)
InB <- sapply(allttimest2,PHealthyToCom1and2)
InDeath <- 1- InC1 - InC2 - InH - InB

if(SavePlot){pdf("Prop-in-each-state.pdf",width=9,height=7)}

par(mfrow=c(1,1),mai=c(0.9,0.9,0.4,0.4))
plot(allttimest2,InDeath,type="l",ylim=c(0,1.1),axes=FALSE,xlab="time t (years)",ylab="Probability of being in the state at t")
polygon(c(allttimest2,rev(allttimest2)),c(rep(0,length(allttimest2)),rev(InDeath)),col=1)

lines(allttimest2,InDeath+InC1,type="l",col=col1)
polygon(c(allttimest2,rev(allttimest2)),c(InDeath,rev(InDeath+InC1)),col=col1)

lines(allttimest2,InDeath+InC1+InC2,type="l",col=col2)
polygon(c(allttimest2,rev(allttimest2)),c(InDeath+InC1,rev(InDeath+InC1+InC2)),col=col2)

lines(allttimest2,InDeath+InC1+InC2+InB,type="l",col=col12)
polygon(c(allttimest2,rev(allttimest2)),c(InDeath+InC1+InC2,rev(InDeath+InC1+InC2+InB)),col=col12)

lines(allttimest2,InDeath+InC1+InC2+InB+InH,type="l",col=col0)
polygon(c(allttimest2,rev(allttimest2)),c(InDeath+InC1+InC2+InB,rev(InDeath+InC1+InC2+InB+InH)),col=col0)

legend("top",fill=c("black",col1,col2,col12,col0),legend=c("Death","Comorb 1","Comorb 2","Comorb 1 & 2","Healthy"),ncol=5)
axis(1)
axis(2,las=2,at=c(0,0.25,0.5,0.75,1),labels=c(0,25,50,75,100))


if(SavePlot){dev.off()}
# }}}


#----------------------------------------------------------------------
### HRandRR.R ends here
