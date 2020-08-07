### Generate-data.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Jun 29 2020 (14:02) 
## Version: 
## Last-Updated: Jul  1 2020 (16:13) 
##           By: Paul Blanche
##     Update #: 46
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

GenerateData <- function(
                         #--- hazards 
                         a0=0.15/(5/1.25),  # baselline (i.e. hazard in healthy group)
                         a01=1/(5/1.25),    # hazard for moving from healthy to comorbidity 1
                         a02=0.25/(5/1.25), # hazard for moving from healthy to comorbidity 2
                         a12=0.25/(5/1.25), # hazard for moving from comorbidity 1 to comorb 1 and 2 
                         a21=0.01/(5/1.25), # hazard for moving from comorbidity 2 to comorb 2 and 1 
                         HR1=6, # Hazard ratio of the time dep Cox model, to model the increase in hazard due to comorbidity 1
                         HR2=2, # Hazard ratio of the time dep Cox model, to model the increase in hazard due to comorbidity 2
                         #--- other parameters
                         n=1000,          # sample size
                         myseed=20191121, # seed 
                         ToScaleTime=5/1.25  # trick to make the time unit nicer
                         ){
    set.seed(myseed)
    # First generate exponential (latent) times.
    # Note: this is just a convenient programming trick, of course we do not assume such times exist.
    T0death <- rexp(n,rate=a0)
    T01 <- rexp(n,rate=a01)
    T02 <- rexp(n,rate=a02)
    T12 <- rexp(n,rate=a12)
    T1death <- rexp(n,rate=a0*HR1)
    T21 <- rexp(n,rate=a21)
    T2death <- rexp(n,rate=a0*HR2)
    T12death <- T21death <-rexp(n,rate=a0*HR2*HR1)
    mat.all.times <- cbind(T0death,
                           T01,
                           T02,
                           T12,
                           T1death,
                           T21,
                           T2death,
                           T12death)
    ## Define life trajectories
    trajectory <- time.death <- time.com1 <- time.com2 <- rep(NA,n)
    ## From healthy state
    DiedDirectly <- which((T0death < T01) & (T0death< T02))
    FirstComob1 <- which((T01 < T0death) & (T01 < T02))
    FirstComob2 <- which((T02 < T0death) & (T02 < T01)) 
    ## From comob 1 state
    DiedComob1 <- which( ((T01 < T0death) & (T01 < T02)) & (T1death < T12)) 
    Comob2AfterComob1 <- which( ((T01 < T0death) & (T01 < T02)) & (T12 < T1death)) 
    ## From comob 2 state
    DiedComob2 <- which( ((T02 < T0death) & (T02 < T01)) & (T2death < T21)) 
    Comob1AfterComob2 <- which( ((T02 < T0death) & (T02 < T01)) & (T21 < T2death)) 
    ## From comob 1 and 2
    DiedComob1and2 <- c(Comob1AfterComob2,Comob2AfterComob1)
    ##
    ## ------ define time to events-----
    ##
    trajectory[DiedDirectly] <- "Died Directly"
    time.death[DiedDirectly] <- T0death[DiedDirectly]
    time.com1[DiedDirectly] <- time.com2[DiedDirectly] <- NA
    ##--
    trajectory[DiedComob1] <- "Died with Comorb 1 only"
    time.death[DiedComob1] <- T01[DiedComob1] + T1death[DiedComob1]
    time.com1[DiedComob1] <- T01[DiedComob1]
    time.com2[DiedComob1] <- NA
    ## --
    trajectory[DiedComob2] <- "Died with Comorb 2 only"
    time.death[DiedComob2] <- T02[DiedComob2] + T2death[DiedComob2]
    time.com1[DiedComob2] <- NA
    time.com2[DiedComob2] <- T02[DiedComob2]
    ## --
    trajectory[intersect(FirstComob1,DiedComob1and2)] <- "Died with Comorb 1 (first) and (then) 2"
    time.death[intersect(FirstComob1,DiedComob1and2)] <- T01[intersect(FirstComob1,DiedComob1and2)] + T12[intersect(FirstComob1,DiedComob1and2)] + T12death[intersect(FirstComob1,DiedComob1and2)]
    time.com1[intersect(FirstComob1,DiedComob1and2)] <- T01[intersect(FirstComob1,DiedComob1and2)]
    time.com2[intersect(FirstComob1,DiedComob1and2)] <- T01[intersect(FirstComob1,DiedComob1and2)] + T12[intersect(FirstComob1,DiedComob1and2)]
    ## --
    trajectory[intersect(FirstComob2,DiedComob1and2)] <- "Died with Comorb 2 (first) and (then) 1"
    time.death[intersect(FirstComob2,DiedComob1and2)] <- T02[intersect(FirstComob2,DiedComob1and2)] + T21[intersect(FirstComob2,DiedComob1and2)] + T21death[intersect(FirstComob2,DiedComob1and2)]
    time.com1[intersect(FirstComob2,DiedComob1and2)] <- T02[intersect(FirstComob2,DiedComob1and2)]+T21[intersect(FirstComob2,DiedComob1and2)]
    time.com2[intersect(FirstComob2,DiedComob1and2)] <- T02[intersect(FirstComob2,DiedComob1and2)]
    #-- data to return ----
    cbind.data.frame(id=1:n,
                     time.death,
                     time.com1,
                     time.com2,
                     trajectory
                     )   
}
message("\nGenerate-data.R: loaded.\n")
#----------------------------------------------------------------------
### Generate-data.R ends here
