## Compute Ferret TH maps (dopaminergic amacrines)
## ---- fth-packages
## To run standalone we need this initial chunk
library(splancs)

## ---- fth-init
## This version requires splancs.

kl <- function(k, show.l=TRUE) {
  ## Return K or L function for plotting.
  if (show.l)
    sqrt(k/pi)
  else
    k
}

myKenv.labelk <- function(pts1, pts2, poly, nsim, s, quiet=FALSE, plot=FALSE) {
  ## This function based on Kenv.label from splancs library.

  
  score.var <- function(k1, k2, k12) {
    ## Score taken from (Diggle, 1986, bottom p122)
    sum(apply( rbind( sqrt(k1), sqrt(k2), sqrt(k12)), 2, var))
  }

  k1.hat <- khat(pts1, poly, s)
  k2.hat <- khat(pts2, poly, s)
  k12.hat<- k12hat(pts1, pts2, poly, s)
  u <- rep(0, length=nsim+1)
  
  kmax <- rep(0, length = length(s))
  kmin <- rep(1e+34, length = length(s))

  u[1] <- score.var(k1.hat, k2.hat, k12.hat)
  
  for (isim in (1:nsim)) {
    if (!quiet) 
      cat("Doing labelling ", isim, "/", nsim, "\n")
    labpts <- rLabel(pts1, pts2)
    k1sim <- khat(labpts[[1]], poly, s)
    k2sim <- khat(labpts[[2]], poly, s)
    k12sim <- k12hat(labpts[[1]], labpts[[2]], poly, s)
    u[isim+1] <- score.var(k1sim, k2sim, k12sim)    
    
    diffk <- kl(k12sim)
    kmax <- pmax(kmax, diffk)
    kmin <- pmin(kmin, diffk)
  }

  if (plot) {
    plot(s, kl(k12.hat), type="l", col="red", lty=2,
         xlab=expression(paste("distance (", mu, "m)")), 
         ylab=expression(L[12]),
         ylim = range(s))
    lines(s, kmin); lines(s, kmax)
  }
  
  list(lower = kmin, upper = kmax, real=k12.hat, u=u, rank=rank(u)[1])
}


######################################################################
## End of functions
######################################################################


## TODO: get the data from the package, not from my homedir!!
f9942.1 <- matrix(scan("~/mosaics/ferret-th/data/f9942i.txt"),
                  ncol=2, byrow=T)
f9942.2 <- matrix(scan("~/mosaics/ferret-th/data/f9942g.txt"),
                  ncol=2, byrow=T)

fth.1 <- f9942.1
fth.2 <- f9942.2
fth.w <- cbind(0, 2500, 0, 2500)
fth.steps <- seq(from=0, to=2500/4, length=150)

## ---- plot-bivariate-fthmaps

plot.fth.maps = function(fth.1, fth.2, relabel=F) {
  fth.soma.rad = 22

  if (relabel) {
    n1 = nrow(fth.1)
    n2 = nrow(fth.2)
    fth.0 = rbind(fth.1, fth.2)
    new.type.1 = sample(n1+n2, n1, replace=FALSE)
    fth.1 = fth.0[new.type.1,]
    fth.2 = fth.0[-new.type.1,]
  }
  
  plot(fth.1, asp=1, type='n', xlab='', ylab='',
       bty='n', xaxt='n', yaxt='n')
  symbols(fth.1[,1], fth.1[,2], circles=rep(fth.soma.rad, nrow(fth.1)),
          lwd=0.2,
          inches=FALSE, add=TRUE)
  symbols(fth.2[,1], fth.2[,2], circles=rep(fth.soma.rad, nrow(fth.2)),
          inches=FALSE, add=TRUE,bg='black')
  rect(fth.w[1], fth.w[3], fth.w[2], fth.w[4])
}

pdf(file='fthmaps.pdf', width=6, height=4)
real.field = 5
par(mar=c(1,1,1,1))
par(mfrow=c(2,3))
for(i in 1:6) {
  plot.fth.maps(fth.1, fth.2, relabel= (i!=real.field))
  title(main=i)
}
dev.off()


## ---- fth-k12
## plot the K function
pts1 = f9942.1
pts2 = f9942.2
steps = fth.steps
datapoly = bboxx(bbox(matrix(fth.w,2,2)))
nsims = 99
par(las=1, bty='n')
klab <- myKenv.labelk(pts1, pts2, datapoly, nsims, steps,
                      quiet=T, plot=TRUE)
rank = klab$rank
pval = rank/(nsims+1)
title(main=sprintf("rank %d p %.2f", rank, pval))

