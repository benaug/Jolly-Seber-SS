NimModel <- nimbleCode({
  lambda.g1 ~ dunif(0,4000) #Expected starting population size
  N[1] ~ dpois(lambda.g1) #Realized starting population size
  U[1] <- N[1] #number of unmarked inds in pop on occasion 1
  M[1] <- 0 #number of marked inds in pop on occasion 1
  
  ##Recruitment and Survival (with removals)##
  for(g in 1:(G-1)){ #loop over periods
    M.plus[g] <- M[g] - m[g] + R[g] #number of marked individuals added after capture/removal
    
    #per capita recruitment (recruited as unmarked)
    gamma[g] ~ dunif(0,2)
    EB[g] <- N[g]*gamma[g] #expected births/recruits
    B[g] ~ dpois(EB[g]) #births/recruits
    
    #survival for unmarked
    phi[g] ~ dunif(0,1)
    #Schofield used this, but I think it is wrong when you remove unmarked inds on first capture
    # U.prime[g] ~ dbinom(size=N[g] - M.plus[g],p=phi[g]) #surviving unmarked from g to g+1
    #this is correct.
    U.prime[g] ~ dbinom(size=U[g] - u[g],p=phi[g]) #surviving unmarked from g to g+1
    
    #survival for marked - factor into known and unknown survival events
    #1) known survival events - prop to phi^T, use ones trick
    T.prob[g] <- pow(phi[g],T[g])
    dummy.data.T[g] ~ dbern(T.prob[g])
    
    #2) unknown survival events - Number in pop + number released (now marked) - known survivals
    S[g] ~ dbinom(size= M.plus[g] - T[g],p=phi[g]) #number survive between g and g+1
    
    #N, U, and M in next period
    U[g+1] <- U.prime[g] + B[g] #unmarked in g+1 is unmarked that survived g plus births from g
    M[g+1] <- S[g] + T[g] #number marked in next session is unknown plus known survivors
    N[g+1] <- U[g+1] + M[g+1]
  }
  
  #Detection
  for(g in 1:G){ #unmarked individuals
    p[g] ~ dunif(0,1)
    u[g] ~ dbinom(p=p[g],size=U[g])
  }
  for(g in 2:G){ #marked individuals
    # Schofield used ones trick. underflow if too many detections
    # m.prob[g] <- pow(p[g],m[g])*pow((1-p[g]),M[g]-m[g])
    # dummy.data.m[g] ~ dbern(m.prob[g])
    m[g] ~ dMarkedDetect(p=p[g],M=M[g]) #distribution that avoids underflow
  }
})

#distribution for detection of marked individuals
dMarkedDetect <- nimbleFunction(
  run = function(x = double(0), p = double(0), M = double(0),log = integer(0)) {
    returnType(double(0))
    lp <- x*log(p) + (M-x)*log(1-p)
    return(lp)
  }
)

#make random data generator to make nimble happy
rMarkedDetect <- nimbleFunction(
  run = function(n = integer(0),p = double(0), M = double(0)) {
    returnType(double(0))
    out <- rbinom(1,size=M,p=p)
    return(out)
  }
)