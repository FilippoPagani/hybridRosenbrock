###########################
#### Hybrid Rosenbrock ####
###########################




## set parameters, challenging case ##
a=1/20 
b=c=100/20
mu=1

## set parameters, easier case ##
#a=50/20
#b=c=100/20
#mu=.2




#### Visualise Distribution ####

n=1e+4 # n. samples

# block 1 #
x=rnorm(n, mu , 1/sqrt(2*a)) # x is common to all blocks
y=rnorm(n, x^2, 1/sqrt(2*b) )
z=rnorm(n, y^2, 1/sqrt(2*c) )

# block 2 (including x, sampled above) #
u=rnorm(n, x^2, 1/sqrt(2*b) )
v=rnorm(n, u^2, 1/sqrt(2*c) )

resDir=matrix(c(x,y,z,u,v),length(x),5)
rm(x,y,z,u,v)

pairs(resDir, upper.panel = NULL) #standard R function





#### sample with Random Walk ####

rw <- function(target, N, x, step, Sigma=diag(length(x)), thin=1)
{
  ptm=proc.time()
  #alpha=1e+6
  cat("Tot run = N*thin = ",N*thin,"\n")
  if ( sum(Sigma - diag(length(x)) )==0  ){
    ident=TRUE
  } else {
    ident=FALSE
    Sigma.chol <- t( chol(Sigma) )
  }
  dd=length(x) 
  acc=0
  samples = matrix(0, N, length(x) )
  samples[1,] <- x
  run=N*thin
  target_prop = target_x = target(x)
  for (i in 2:run )
  {
    if(i%%(run/10)==0) cat("Progress:",i/run*100,"% \n") #track progress
    
    if (ident) prop <- x + step*rnorm(dd)
    else prop = x + step * Sigma.chol %*% rnorm(dd)
    
    target_prop = target(prop)
    
    if (runif(1) < exp( target_prop - target_x ) ){
      x <- prop
      acc = acc + 1
      target_x = target_prop
    }
    
    if( (i%%thin)==0 ) samples[i/thin,] = x
  }
  time=proc.time()-ptm
  cat('Run Time',round(time[1]/60,digits=2),'Min','\n')
  #cat("Acc Rate",acc/run,"\n")
  system(paste("echo Acc",acc/run))
  return(list(samples,time[1]))
}

logTarget6=function(q){-a*(q[1]-mu)^2 -b*(q[2] - q[1]^2)^2 -b*(q[3] - q[2]^2)^2 -b*(q[4] - q[1]^2)^2 -b*(q[5] - q[4]^2)^2 } #U=+log
resRW=rw(logTarget6, 100000, matrix(1,5,1), .2)




#### compare RWM results with iid results ####

pairs(resRW[[1]], upper.panel = NULL)





#### END ####