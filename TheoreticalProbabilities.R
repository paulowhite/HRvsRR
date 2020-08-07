### TheoreticalProbabilities.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Jul  1 2020 (15:24) 
## Version: 
## Last-Updated: Jul  2 2020 (10:07) 
##           By: Paul Blanche
##     Update #: 89
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# Probability to be alive with comorbidity 1 only in t years, given being currently alive with no comorbidity.
PHealthyToCom1 <- function(t,
                           a0=0.15/(5/1.25),
                           a01=1/(5/1.25),
                           a02=0.25/(5/1.25),
                           a12=0.25/(5/1.25),
                           HR1=6){
                               ((a01*exp(-(a0*HR1 + a12)*t))/(a01+a02+a0-a12-a0*HR1))*(1-exp(-(a01+a02+a0-a12-a0*HR1)*t))
}
# Probability to be alive with comorbidity 1 only in t years, given being currently alive with no comorbidity.
PHealthyToCom2 <- function(t,
                           a0=0.15/(5/1.25),  
                           a01=1/(5/1.25),    
                           a02=0.25/(5/1.25), 
                           a21=0.01/(5/1.25), 
                           HR2=2){
                               ((a02*exp(-(a0*HR2 + a21)*t))/(a02+a01+a0-a21-a0*HR2))*(1-exp(-(a02+a01+a0-a21-a0*HR2)*t))
}
# Probability to remain alive without comorbidity within t years.
PStillHealthy <- function(t,
                          a0=0.15/(5/1.25),  
                          a01=1/(5/1.25),    
                          a02=0.25/(5/1.25)
                          ){
    exp(-(a0+a01+a02)*t)
}
# Risk of dying within t years given being alive with comorbidity 1 only.
Risk1 <- function(t,
                  a0=0.15/(5/1.25),  
                  a12=0.25/(5/1.25), 
                  HR1=6, 
                  HR2=2){
    1- exp(-t*(a0*HR1+a12))-( a12*exp(-a0*HR1*HR2*t)/(a0*HR1 + a12-a0*HR1*HR2))*(1-exp(-t*(a0*HR1+a12-a0*HR1*HR2)))
}
# Probability of being alive with both comorbidities in t years, given being alive with comorbidity 1 only.
PCom1toCom1and2 <- function(t,
                            a0=0.15/(5/1.25),  
                            a12=0.25/(5/1.25), 
                            HR1=6, 
                            HR2=2){
                                ((a12*exp(-a0*HR1*HR2*(t)))/(a0*HR1 + a12 - a0*HR1*HR2))*( 1 - exp(-(a0*HR1 + a12 - a0*HR1*HR2)*(t))  )
}
# Probability of still being alive with comorbidity 1 only in t years.
PStillOnlyCom1 <- function(t,
                           a0=0.15/(5/1.25),  
                           a12=0.25/(5/1.25), 
                           HR1=6){
    exp(-(a0*HR1 + a12)*t)
}
# Risk of dying within t years given being alive with comorbidity 2 only.
Risk2 <- function(t,
                  a0=0.15/(5/1.25),  
                  a21=0.01/(5/1.25), 
                  HR1=6, 
                  HR2=2){
    1- exp(-t*(a0*HR2+a21)) - ( a21*exp(-a0*HR1*HR2*t)/(a0*HR2 + a21-a0*HR1*HR2))*(1-exp(-t*(a0*HR2+a21-a0*HR1*HR2)))
}
# Probability of being alive with both comorbidities in t years, given being alive with comorbidity 2 only.
PCom2toCom1and2 <- function(t,
                            a0=0.15/(5/1.25),  
                            a21=0.01/(5/1.25), 
                            HR1=6, 
                            HR2=2){
                                ((a21*exp(-a0*HR1*HR2*(t)))/(a0*HR2 + a21 - a0*HR1*HR2))*( 1 - exp(-(a0*HR2 + a21 - a0*HR1*HR2)*t) )
}
# Probability of still being alive with comorbidity 2 only in t years.
PStillOnlyCom2 <- function(t,
                           a0=0.15/(5/1.25),  
                           a21=0.01/(5/1.25), 
                           HR2=2){
    exp(-(a0*HR2 + a21)*t)
}
# Probability to be alive with both comorbidities in t years, given being currently alive with none.
PHealthyToCom1and2 <- function(t,
                               a0=0.15/(5/1.25),  
                               a01=1/(5/1.25),    
                               a02=0.25/(5/1.25), 
                               a12=0.25/(5/1.25), 
                               a21=0.01/(5/1.25), 
                               HR1=6, 
                               HR2=2){
    term12.1 <- (a12*a01*exp(-a0*HR1*HR2*t))/(a0*HR1 + a12 - a0*HR1*HR2)
    term12.2 <- (1-exp(-(a0+a01+a02-a0*HR1*HR2)*t))/(a0 + a01 + a02 - a0*HR1*HR2)
    term12.3 <- exp(-(a0*HR1+a12-a0*HR1*HR2)*t)*(1-exp(-(a0+a01+a02-a0*HR1-a12)*t))/(a0+a01+a02-a0*HR1-a12)
    Delta12 <- term12.1*(term12.2 - term12.3)
    ##
    term21.1 <- (a21*a02*exp(-a0*HR1*HR2*t))/(a0*HR2 + a21 - a0*HR1*HR2)    
    term21.2 <- term12.2    
    term21.3 <- exp(-(a0*HR2+a21-a0*HR1*HR2)*t)*(1-exp(-(a0+a01+a02-a0*HR2-a21)*t))/(a0+a01+a02-a0*HR2-a21)    
    Delta21 <- term21.1*(term21.2 - term21.3)
    ##
    Delta12+Delta21   
}
# Risk of dying within t years given being alive with both comorbidities.
Risk1and2 <- function(t,
                      a0=0.15/(5/1.25),  
                      HR1=6, 
                      HR2=2){
    1- exp(-(a0*HR2*HR1)*t)
}

#----------------------------------------------------------------------
### TheoreticalProbabilities.R ends here
