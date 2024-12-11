#this is an example of how to compute the sufficient statistics from the capture history and number removed per occasion.

source("sim.JS.removal.R")

G <- 20 #number of periods
lambda.g1 <- 250 #expected N in year 1
gamma <- rep(0.2,G-1) #yearly per-capita recruitment
phi <- rep(0.75,G-1)
p <- rep(0.6,G) #yearly detection probabilities
#parameters to simulate removals, but not required for model estimation
p.remove.marked <- rep(0,G) #probability captured unmarked individuals are removed
p.remove.unmarked <- rep(0.2,G) #probability captured marked individuals are removed

data <- sim.JS.removal(lambda.g1=lambda.g1,gamma=gamma,
                       phi=phi,p=p,p.remove.marked=p.remove.marked,
                       p.remove.unmarked=p.remove.unmarked,G=G)

#compute sufficient statistics from capture history and history for removals
y.obs <- data$y.obs
N.removed <- data$N.removed

#compute marked and unmarked number of captures, u and m
n.cap <- nrow(y.obs)
y.unmarked <- array(0,dim=c(n.cap,G))
for(i in 1:n.cap){
  position <- Position(y.obs[i,],f=function(x){x>0})#occasion of first capture, NA if no capture
  y.unmarked[i,position] <- y.obs[i,position]
}
y.marked <- y.obs
y.marked[y.unmarked>0] <- 0

u <- colSums(y.unmarked>0)
m <- colSums(y.marked>0)

all(u==data$u)
all(m==data$m)

#compute releases, R
R <- m + u #number caught in sample g and released (ignoring removals)
R <- R - N.removed #subtract off removals
R <- R[-G] #remove last year

all(R==data$R)

#compute T, number of known survival events from marked individuals
z.known <- y.obs*0
for(i in 1:nrow(y.obs)){
  idx  <- which(y.obs[i,]==1) #nothing to fill in
  if(length(idx)>1){
    z.known[i,(min(idx)+1):max(idx)] <- 1 #ad one because you could be recruit
  }else{#if you're a singleton, no z fixed, could be recruit
    z.known[i,] <- 0
  }
}

T <- colSums(z.known)[-1] #known survivals of marked inds

all(T==data$T)

#sufficient statistics:
u
m
R
T


