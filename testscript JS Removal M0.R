#this model/script handles the cases of no removals, unmarked removals, marked removals, or both
#M_0: all parameters are fixed over occasions. See M_t testscript for occasion effects

library(nimble)
library(coda)

source("sim.JS.removal.R")
source("Nimble Model JS Removal M0.R")

G <- 20 #number of periods
lambda.g1 <- 250 #expected N in year 1
gamma <- rep(0.2,G-1) #yearly per-capita recruitment
phi <- rep(0.75,G-1)
p <- rep(0.6,G) #yearly detection probabilities
#parameters to simulate removals, but not required for model estimation
p.remove.marked <- rep(0,G) #probability captured unmarked individuals are removed
p.remove.unmarked <- rep(0,G) #probability captured marked individuals are removed

data <- sim.JS.removal(lambda.g1=lambda.g1,gamma=gamma,
               phi=phi,p=p,p.remove.marked=p.remove.marked,
               p.remove.unmarked=p.remove.unmarked,G=G)

data$truth$N
data$truth$M.removed
data$truth$U.removed

#init from truth (for testing/development)
# N.init <- data$truth$N
# D.init <- data$truth$D
# M.init <- data$truth$M
# U.init <- data$truth$U
# B.init <- data$truth$B
# S.init <- data$truth$S
# M.plus.init <- data$truth$M.plus
# U.prime.init <- data$truth$U.prime

# init from data (may not be the most efficient way, but it works)
M.init <- rep(NA,G)
D.init <- rep(NA,G-1)
M.init[1] <- 0
D.init[1] <- M.init[1]-data$m[1]+data$R[1]-data$T[1]
for(g in 2:G){
  M.init[g] <- M.init[g-1] - data$m[g-1] + data$R[g-1] - D.init[g-1]
  D.init[g] <-  (M.init[g] - data$m[g] + data$R[g]) - data$T[g] #inds that could die minus known survivals
}
D.init <- D.init[-G]

U.init <- data$u + 10
N.init <- M.init + U.init
U.prime.init <- U.init[-G] - data$u[-G]
B.init <- U.init[-1] - U.prime.init
M.plus.init <- M.init[-G] - data$m[-G] + data$R
S.init <- M.plus.init-data$T

#tests
all(M.init[2:G] == (M.init[1:(G-1)] - data$m[1:(G-1)] + data$R[1:(G-1)] - D.init[1:(G-1)]))
all(M.init[-1] >= data$m[-1])
all(data$u<=U.init)
all(U.prime.init<=(U.init[-G]-data$u[-G]))
all((M.plus.init-data$T)>=S.init)

#must be equal
all(N.init[-1]==(U.init[-1] + M.init[-1]))
all(U.init[-1] == (U.prime.init + B.init))
all(M.init[-1] == M.init[-G] - data$m[-G] + data$R - D.init)

#get p, phi, gamma inits from initialized latent variables
#if you have a lot of data, you may need to use these for starting logprobs to be finite
#might as well always use them
p.init <- (data$m+data$u)/(M.init+U.init)
phi.init <- (S.init+data$T+U.prime.init)/N.init[-G]
gamma.init <- B.init/N.init[-G]

#use yearly inits above for yearly parameters, average them for fixed
p.init <- mean(p.init)
phi.init <- mean(phi.init)
gamma.init <- mean(gamma.init)

constants <- list(G=G) 
Niminits <- list(p=p.init,gamma=gamma.init,phi=phi.init,lambda.g1=N.init[1],
                 N=N.init,S=S.init,M=M.init,U=U.init,B=B.init,U.prime=U.prime.init)
Nimdata <- list(u=data$u,m=data$m,R=data$R,T=data$T,dummy.data.T=rep(1,G-1))

# set parameters to monitor
parameters <- c('lambda.g1','gamma','p','phi')
nt <- 10 #thinning rate
parameters2 <- c('N')
nt2 <- 10

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,monitors2=parameters2, thin2=nt2)
#probably don't need block sampler, but can use one
# conf$removeSampler('gamma','p','phi')
# conf$addSampler(target = c("phi","gamma","p"),type = 'AF_slice',
#                 control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmodel$calculate()

#if starting logprob not finite, fix inits. check these below to find problem init.
# Cmodel$logProb_B
# Cmodel$logProb_U.prime
# Cmodel$logProb_S
# Cmodel$logProb_dummy.data.T
# Cmodel$logProb_m  #if too many captures, can get underflow here. need to fix this.
# Cmodel$logProb_u

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(10000,reset=TRUE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

burnin <- 100
plot(mcmc(mvSamples[-c(1:burnin),]))

round(cor(mcmc(mvSamples[-c(1:burnin),])),2)


mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
burnin <- 100
plot(mcmc(mvSamples2[-c(1:burnin),]))

data$truth$N

