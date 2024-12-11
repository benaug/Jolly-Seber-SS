sim.JS.removal <- function(lambda.g1=NA,gamma=NA,phi=NA,p=NA,
                           p.remove.marked=NA,p.remove.unmarked=NA,G=NA){
  N <- M <- U  <- rep(NA,G)
  N.recruit <- N.survive <- N.survive.M <- N.survive.U <- D <- ER <- U.prime <- rep(NA,G-1)
  
  N[1] <- rpois(1,lambda.g1)
  M[1] <- 0
  U[1] <- N[1]
  #Easiest to increase dimension of z as we simulate bc size not known in advance.
  #with removals, need to simulate popdy and capture iteratively
  z <- y <- y.U.removed <- y.M.removed <- matrix(0,N[1],G)
  status <- matrix(0,N[1],G) #start in unmarked (1), transition to marked (2) after capture
  status[,1] <- 1 #start as unmarked
  z[1:N[1],1] <- 1
  N.super <- N[1]
  
  #occasion 1
  #capture
  y[,1] <- rbinom(nrow(z),1,p[1]*z[,1])
  #removal
  y.U.removed[,1] <- rbinom(N.super,1,y[,1]*p.remove.unmarked[1]) #all unmarked
  
  #between occasion 1 and 2
  #recruits
  ER[1] <- N[1]*gamma[1]
  N.recruit[1] <- rpois(1,ER[1])
  #increase dimensions
  status[,2] <- status[,1] #shift current status to next session for alive guys, update below
  if(N.recruit[1]>0){
    z <- rbind(z,matrix(0,nrow=N.recruit[1],ncol=G))
    z[(N.super+1):(N.super+N.recruit[1]),2] <- 1 #enter as alive, z=1
    status <- rbind(status,matrix(0,nrow=N.recruit[1],ncol=G))
    status[(N.super+1):(N.super+N.recruit[1]),2] <- 1 #enter as unmarked=1
    y.U.removed <- rbind(y.U.removed,matrix(0,nrow=N.recruit[1],ncol=G))
    y.M.removed <- rbind(y.M.removed,matrix(0,nrow=N.recruit[1],ncol=G))
    y <- rbind(y,matrix(0,nrow=N.recruit[1],ncol=G))
  }
  
  #survival
  idx <- which(z[,1]==1&y.U.removed[,1]==0&y.M.removed[,1]==0) #remove removals from survival. recruits not used bc their z[,1]=0
  z[idx,2] <- rbinom(length(idx),1,phi[1])
  U.prime[1] <- sum(z[,1]==1&z[,2]==1&status[,1]==1&y[,1]==0) #unmarked survivors, cant remove if not captured
  N.survive[1] <- sum(z[,1]==1&z[,2]==1)
  N[2] <- N.recruit[1] + N.survive[1]
  N.super <- N.super + N.recruit[1]
  
  #update status
  status[which(y[,1]==1),2] <- 2 #if captured, marked in next period
  status[which(y.U.removed[,1]==1),2] <- 0 #removed
  status[which(y.M.removed[,1]==1),2] <- 0 #removed
  N.survive.M[1] <-  sum(z[,1]==1&z[,2]==1&status[,2]==2)
  N.survive.U[1] <-  sum(z[,1]==1&z[,2]==1&status[,2]==1)
  
  D[1] <- sum(z[,1]==1&z[,2]==0&status[,2]==2) #number marked deaths
  U[2] <- sum(status[,2]==1&z[,2]==1)
  M[2] <- sum(status[,2]==2&z[,2]==1)
  
  for(g in 2:(G-1)){
    #on occasion g
    #capture
    y[,g] <- rbinom(nrow(z),1,p[g]*z[,g])
    #removal
    y.U.removed[,g] <- rbinom(N.super,1,p.remove.unmarked[g]*(y[,g]==1&status[,g]==1))
    y.M.removed[,g] <- rbinom(N.super,1,p.remove.marked[g]*(y[,g]==1&status[,g]==2))
    
    #between occasions g and g+1
    
    #Simulate recruits
    ER[g] <- N[g]*gamma[g]
    N.recruit[g] <- rpois(1,ER[g])
    #add recruits to z, y, y.removed
    status[,g+1] <- status[,g] #shift current status to next session then update
    if(N.recruit[g]>0){
      z <- rbind(z,matrix(0,nrow=N.recruit[g],ncol=G))
      z[(N.super+1):(N.super+N.recruit[g]),g+1] <- 1
      y.M.removed <- rbind(y.M.removed,matrix(0,nrow=N.recruit[g],ncol=G))
      y.U.removed <- rbind(y.U.removed,matrix(0,nrow=N.recruit[g],ncol=G))
      y <- rbind(y,matrix(0,nrow=N.recruit[g],ncol=G))
      status <- rbind(status,matrix(0,nrow=N.recruit[g],ncol=G))
      status[(N.super+1):(N.super+N.recruit[g]),g+1] <- 1 #enter as unmarked=1
    }
    
    #survival
    idx <- which(z[,g]==1&y.U.removed[,g]==0&y.M.removed[,g]==0) #remove removals from survival
    z[idx,g+1] <- rbinom(length(idx),1,phi[g])
    U.prime[g] <- sum(z[,g]==1&z[,g+1]==1&status[,g]==1&y[,g]==0) #unmarked survivors 
    N.survive[g] <- sum(z[,g]==1&z[,g+1]==1)
    N[g+1] <- N.recruit[g] + N.survive[g]
    N.super <- N.super + N.recruit[g]
    
    #update status
    status[which(y[,g]==1),g+1] <- 2 #if captured, marked in next period
    status[which(status[,g]==2),g+1] <- 2 #if previously marked, still marked
    status[which(y.U.removed[,g]==1),g+1] <- 0 #removed
    status[which(y.M.removed[,g]==1),g+1] <- 0 #removed
    N.survive.M[g] <-  sum(z[,g]==1&z[,g+1]==1&status[,g+1]==2)
    N.survive.U[g] <-  sum(z[,g]==1&z[,g+1]==1&status[,g+1]==1)
    D[g] <- sum(z[,g]==1&z[,g+1]==0&status[,g+1]==2) #number marked deaths
    U[g+1] <- sum(status[,g+1]==1&z[,g+1]==1)
    M[g+1] <- sum(status[,g+1]==2&z[,g+1]==1)
  }
  if(any(N.recruit+N.survive!=N[2:G]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")
  #final occasion capture
  y[,G] <- rbinom(nrow(z),1,p[g]*z[,G])
  
  #discard individuals never detected
  y.obs <- y[rowSums(y)>0,]
  z.obs <- z[rowSums(y)>0,]
  n.cap <- nrow(y.obs)
  y.unmarked <- array(0,dim=c(n.cap,G))
  for(i in 1:n.cap){
    position <- Position(y.obs[i,],f=function(x){x>0})#occasion of first capture, NA if no capture
    y.unmarked[i,position] <- y.obs[i,position]
  }

  y.marked <- y.obs
  y.marked[y.unmarked>0] <- 0
  U.removed <- colSums(y.U.removed)
  M.removed <- colSums(y.M.removed)
  u <- colSums(y.unmarked>0)
  m <- colSums(y.marked>0)
  R <- m + u #number caught in sample g and released (ignoring removals)
  R <- R - U.removed - M.removed #subtract off removals
  R <- R[-G] #remove last year
  
  #these are marked individuals known to survive. can compute one of two ways
  # T <- rep(NA,G-1) #number inds caught up to and including g that were also captured after g
  # T[1] <- sum(y[,1]>0&rowSums(y[,2:G])>0) #caught first year and in a subsequent year
  # for(g in 2:(G-2)){
  #   T[g] <- sum(rowSums(y[,1:g])>0&rowSums(y[,(g+1):G])>0) #caught in year g or before and in a subsequent year
  # }
  # T[G-1] <- sum(rowSums(y[,1:(G-1)])>0&y[,G]>0) #caught in next to last year or before and in last year
  
  #identify known and unknown survival events for marked inds
  z.known <- y.obs*0
  # S <- rep(0,G-1)
  S <- rep(0,G)
  for(i in 1:nrow(y.obs)){
    idx  <- which(y.obs[i,]==1) #nothing to fill in
    if(length(idx)>1){
      z.known[i,(min(idx)+1):max(idx)] <- 1 #ad one because you could be recruit
    }else{#if you're a singleton, no z fixed, could be recruit
      z.known[i,] <- 0
    }
    if(max(idx)<G){
      S <- S + 1*((1:(G))>max(idx))*(z.obs[i,]==1)
    }
  }
  T <- colSums(z.known)[-1] #known survivals of marked inds
  S <- S[-1]
  # (S+T)==M[-1]
  
  M.plus <- M[-G] - m[-G] + R
  # S <- M.plus - T #unknown survival events of marked inds
  
  B <- N.recruit
  
  # #make sure R is consistent with M and D
  M.curr <- 0
  for(g in 1:(G-1)){
    M.curr <-  M[g] - m[g] + R[g] - D[g]
    if(M.curr!=M[g+1])stop("bug in code")
  }
  
  #p is right
  mean((u+m)/(U+M))
  mean(m[-1]/M[-1])
  mean(u/U)

  #unmarked phi is right
  mean(U.prime/(N[-G]-M.plus[-G]))
  mean(U.prime/(U[-G]-u[-G]))

  #marked phi is right
  #complement of number of marked deaths / opportunities for marked deaths
  1-(sum(D)/sum(M[-G] + (R - m[-G])))
  1-(sum(D)/sum(M[-G] + (R - m[-G]-T)))

  #recruitment is right
  mean(B/N[-G])
  
  # y.removed <- y.U.removed[rowSums(y)>0,] + y.M.removed[rowSums(y)>0,]
  

  #store true data for model building/debugging
  truth <- list(y=y,N=N,N.recruit=N.recruit,N.survive=N.survive,z=z,
                M=M,U=U,U.prime=U.prime,D=D,B=N.recruit,M.plus=M.plus,S=S,
                U.removed=U.removed,M.removed=M.removed)
  

  return(list(u=u,m=m,R=R,T=T,y.obs=y.obs,N.removed=U.removed+M.removed,truth=truth))
}

