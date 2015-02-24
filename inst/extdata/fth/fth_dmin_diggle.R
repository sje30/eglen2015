## Adapted from ~/mosaics/ferret-th/diggle-analysis/dmin_diggle.r

## Run the Diggle analysis on various data sets.  To see which ones
## to run, do (occur "^do") within Emacs.
## Thu 20 Feb 2003
## Can run in Batch using:
## R BATCH dmin_diggle.r

library(splancs)
library(sjedrp)
library(sjedmin)
library(sjevor)

par(bty="n")

kl <- function(k, show.l=TRUE) {
  ## Return K or L function for plotting.
  if (show.l)
    sqrt(k/pi)
  else
    k
}

myKenv.csr <- function (pts, poly, nsim, s, quiet = FALSE, plot = FALSE,
                        my.ylab=1, dmin.mu=NULL, dmin.sd=NULL)  {
  ## Adapted from Kenv.csr from splancs library to do scoring.
  ## If dmin.mu is supplied, we should compare against a dmin simulation
  ## rather than CSR fields.  dmin = 12.3 +- 1.5 for the FTH cells.
  theo <- sqrt(pi * s*s)
  
  score.k <- function(k) {
    ## Return "u" from (Diggle, 1986, p122, top equation).
    sum( (sqrt(k)-theo)^2)
  }
  
  u <- rep(0, length=nsim+1)
  npts <- dim(pts)[1]
  kmax <- rep(0, length = length(s))
  kmin <- rep(1e+34, length = length(s))
  k.real <-  khat(pts, poly, s)
  u[1] <- score.k(k.real)

  ## For generating CSR simulations, do we want to use csr or dmin?
  ## Mon 17 Mar 2003: Checked when CSR constrained by ferret soma size
  ## that results are unaffected; output in fth_dmin_csr1/
  use.dmin <- !is.null(dmin.mu)
  if (use.dmin) 
    rejs <- rep(0, nsim)
  
  for (isim in (1:nsim)) {
    if (!quiet) 
      cat("Doing simulation ", isim, "\n")

    if (use.dmin) {
      bb <- bbox(poly); dwid <- max(bb[,1]) - min(bb[,1])
      dht <- max(bb[,2]) - min(bb[,2])
      dmin <- dminl(dwid, dht, npts, dmin.mu, dmin.sd, quiet=T)
      rejs[isim] <- sum(dmin$nrejects)
      dpoly <- spoints( c(0,0, 0,dht,   dwid,dht,  dwid,0))
      khsim <- khat(as.points( cbind(dmin$x, dmin$y)), dpoly, s)
    } else {
      khsim <- khat(csr(poly, npts), poly, s)
    }
    kmax <- pmax(kmax, khsim)
    kmin <- pmin(kmin, khsim)
    u[isim+1] <- score.k(khsim)
  }

  if (use.dmin) {
    cat("Range of dmin rejects from simulations\n")
    print(summary(rejs))
  }
  
  if(my.ylab == 1) {
    the.ylab <- expression(L[1])
  } else if (my.ylab == 2) {
    the.ylab <- expression(L[2])
  } else {
    the.ylab <- expression(L[1+2])
  }
  
  if(plot) {
    plot(s, kl(k.real), type="l", col="red",
         xlab=expression(paste("distance (", mu, "m)")), 
         ylim=c(0, max(s)),
         ylab=the.ylab)
    lines(s, kl(kmin)); lines(s, kl(kmax))
  }
  
  
  list(lower = kmin, upper = kmax, real=k.real, u=u, rank=rank(u)[1])
}



myKenv.tor <- function (pts1, pts2, poly, nsim, s, quiet = FALSE, plot=FALSE) {
  ## Adapted from Kenv.tor in splancs library to do scoring.
  rect <- bboxx(bbox(poly))
  theo <- sqrt(pi * s*s)
  ## U score taken from top of (Diggle 1986, p122)
  score.k <- function(k) { sum( (sqrt(k)-theo)^2) }
  ##score.k <- function(k) { sum( (1/(s^2)) * (k - theo)^2

                           
  u <- rep(0, length=nsim+1)
  kmax <- rep(0, length = length(s))
  kmin <- rep(1e+34, length = length(s))
  k12.real <-  k12hat(pts1, pts2, poly, s)
  u[1] <- score.k(k12.real)
  for (isim in (1:nsim)) {
    if (!quiet) 
      cat("Doing simulation ", isim, "\n")
    pts2s <- rtor.shift(pts2, poly)
    k12.sim <-  k12hat(pts1, pts2s, poly, s)
    kmax <- pmax(kmax, k12.sim)
    kmin <- pmin(kmin, k12.sim)
    u[isim+1] <- score.k(k12.sim)
  }

  if(plot) {
    plot(s, kl(k12.real), type="l",
         xlab=expression(paste("distance (", mu, "m)")),
         ylab=expression(L[12]),
         col="red", lty=2)
    lines(s, kl(kmin)); lines(s, kl(kmax))
  }
  
  list(lower = kmin, upper = kmax, real=k12.real, u=u, rank=rank(u)[1])
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




diggle.analysis <- function(pts1, pts2, steps, datapoly=NULL,
                            plot=F, controls=F, nsims=199, show.fields=F,
                            do.csr=T, do.tor=T, do.lab=T ) {

  ## Optionally do all the Diggle tests, and return result as list.
  ## Any test result of NULL indicates test not performed.
  ## So, by default, initialise all test results to NULL.
  csr.1 <- csr.2 <- csr.0 <- NULL
  k12.mytor <- klab <- NULL
  
  ## Plot points just to check we have read them in okay.
  if (plot && show.fields) {
    dataset <- paste( deparse(substitute(pts1)), deparse(substitute(pts2)))
    plot(pts1,pch=19,asp=1, main=dataset); points(pts2)
  }

  pts0 <- rbind(pts1,pts2)

  if (is.null(datapoly)) {
    ## No datapoly provided so we should assume it just encloses all the
    ## cells in both fields.
    xmin <- min(pts0[,1]); ymin <- min(pts0[,2])
    xmax <- max(pts0[,1]); ymax <- max(pts0[,2])

    datapoly <- spoints( c(xmin,ymin, xmin, ymax,   xmax, ymax,  xmax,ymin))
  }

  
  ## to remove titles, use this line.
  ##title <- function(x) { print("ignore x");}


  ## Test 1 (Diggle, 1986 page 121), test each component for CSR.

  ## 1a. test type 1 cells against CSR.
  if (do.csr) {
    csr.1 <- myKenv.csr(pts1, datapoly, nsims, steps, quiet=T,plot=plot,
                        my.ylab=1)
    title.csr.1 <- paste("1a. CSR of 1: rank", csr.1$rank,
                         "u", signif(csr.1$u[1],3),
                         "vs [",
                         paste(signif(range(csr.1$u[-1]),3),collapse=","), "]")
    cat(paste("RES: ", title.csr.1, "\n"))
    if (plot)
      title(title.csr.1)
    
    ## 1b. test type 2 cells against CSR
    csr.2 <- myKenv.csr(pts2, datapoly, nsims, steps, quiet=T, plot=T,
                        my.ylab=2)
    title.csr.2 <- paste("1b. CSR of 2: rank", csr.2$rank,
                         "u", signif(csr.2$u[1],3),
                         "vs [",
                         paste(signif(range(csr.2$u[-1]),3),collapse=","), "]")
    cat(paste("RES: ", title.csr.2, "\n"))
    if (plot)
      title(title.csr.2)
    
    
    ## 1c. test type 0 cells against CSR.
    csr.0 <- myKenv.csr(pts0, datapoly, nsims, steps, quiet=T, plot=plot,
                        my.ylab=0)
    title.csr.0 <- paste("1c. CSR of 1+2: rank", csr.0$rank,
                         "u", signif(csr.0$u[1],3),
                         "vs [",
                         paste(signif(range(csr.0$u[-1]),3),collapse=","), "]")
    cat(paste("RES: ", title.csr.0, "\n"))
    if (plot)
      title(title.csr.0)


    ## 1d. CONTROL.  Make N random points, where N is the number of
    ## cells in class 1; these should not be significantly different from
    ## CSR!
    
    if (controls) {
      csr.r <- myKenv.csr(csr(datapoly, dim(pts1)[1]),
                          datapoly, nsims, steps,quiet=T, plot=plot)
      title.csr.r <- paste("1d. Control CSR: rank", csr.r$rank,
                           "u", signif(csr.r$u[1],3), "vs [",
                           paste(signif(range(csr.r$u[-1]),3),collapse=","), "]")
      cat(paste("RES: ", title.csr.r, "\n"))
      if (plot)
        title(title.csr.r)
    }
  }


  if (do.tor) {
    ## Test 2 - toroidal shifting.  Does position of one type influence
    ## position of other?  (Page 122, Diggle 1986).

    ## Test 2a. Test whether 1+2 are spatially independent. 
    k12.mytor <- myKenv.tor(pts1, pts2, datapoly, nsims, steps,
                            quiet=T,plot=plot)
    title.k12tor <- (paste("2a. tor: rank", k12.mytor$rank, "u",
                           signif(k12.mytor$u[1],3), "vs [",
                           paste(signif(range(k12.mytor$u[-1]),3),
                                 collapse=","),"]"))
    
    cat(paste("RES: ", title.k12tor, "\n"))
    if (plot)
      title(title.k12tor)
    
    ## Test 2b. As control, generate two data sets at random, these should
    ## produce a rank less than 95, and be spatially independent.
    if (controls) {
      k12r.mytor <- myKenv.tor(csr(datapoly, dim(pts1)[1]),
                               csr(datapoly, dim(pts2)[1]),
                               datapoly, nsims, steps, quiet=T, plot=plot)
      title.k12r <- paste("2b. Control toroidal - rank",
                          k12r.mytor$rank, "u",
                          signif(k12r.mytor$u[1],3), "vs [",
                          paste(signif(range(k12r.mytor$u[-1]),3),collapse=","),
                          "]")
      cat(paste("RES: ", title.k12r, "\n"))
      if (plot)
        title(title.k12r)
    }
  }


  if (do.lab) {
    ## Test 3 - independence.  Are cell types randomly labelled?
    ## Described middle of page 122, Diggle 1986.
    
    ## Test 3a. See whether type 1 and type 2 are randomly labelled.
    klab <- myKenv.labelk(pts1, pts2, datapoly, nsims, steps,
                         quiet=T, plot=plot)    
    title.randlabel <- paste("3a. rlab: rank", klab$rank, "u",
                             signif(klab$u[1],3), "vs [",
                             paste(signif(range(klab$u[-1]),3),collapse=","),
                             "]" )
    if (plot)
      title(title.randlabel)

    if (!is.null(klab$ri1)) {
      browser()
      stripchart(list(ri1=klab$ri1, ri2=klab$ri2),
                 group.names=c("Class 1", "Class 2"),
                 ylab="regularity index",
                 vertical=T, pch=19, cex=0.5 , method="jitter")
      dx <- 0.1
      segments( 1-dx, klab$ri1[1], 1+dx, klab$ri1[1])
      segments( 2-dx, klab$ri2[1], 2+dx, klab$ri2[1])
    }
    cat(paste("RES: ", title.randlabel, "\n"))

    
    
    ## Test 3b. As control, divide cells into two populations at random.
    if (controls) {
      ron <- sample(dim(pts0)[1], dim(pts1)[1], replace=F)
      rpts1 <- pts0[ron,]; rpts2 <- pts0[-ron,]
      if (plot) {
        plot(rpts1,pch=19,asp=1, main="3b. random labelling of orig data");
        points(rpts2)
      }
      rklab <- myKenv.label(rpts1, rpts2, datapoly, nsims, steps, quiet=T,
                            plot=plot)
      
      title.randlabel.r <- paste("3b. Control rlab: rank",
                                 rklab$rank, "u",
                                 signif(rklab$u[1],3), "vs [",
                                 paste(signif(range(rklab$u[-1]),3),collapse=","),
                                 "]")
      cat(paste("RES: ", title.randlabel.r, "\n"))
      cat(title.randlabel.r)
      if (plot)
        title(title.randlabel.r)
    }
  }

  res <- list(csr.1=csr.1, csr.2=csr.2, csr.0=csr.0,
              tor=k12.mytor, lab=klab, steps=steps)
  res
}

######################################################################
## end of functions
######################################################################

do.ferret.th <- TRUE

if (do.ferret.th) {
  postscript(file="ferret_th_diggle.ps")
  ##par(mfrow=c(2,3), bty="n", las=1)
  par(mfrow=c(2,2), bty="n", las=1)
  f9942.1 <- matrix(scan("~/mosaics/ferret-th/data/f9942i.txt"),
                    ncol=2, byrow=T)
  f9942.2 <- matrix(scan("~/mosaics/ferret-th/data/f9942g.txt"),
                    ncol=2, byrow=T)
  f.steps <- seq(from=0, to=2500/4, length=150)
  f9942.dres <- diggle.analysis(f9942.1, f9942.2,
                                show.fields = TRUE,
                                do.csr=FALSE, plot=T,steps=f.steps)
  
  ## f0102.1 <- matrix(scan("~/mosaics/ferret-th/data/f0102i.txt"),
  ##                   ncol=2, byrow=T)
  ## f0102.2 <- matrix(scan("~/mosaics/ferret-th/data/f0102g.txt"),
  ##                   ncol=2, byrow=T)
  ## f.steps <- seq(from=0, to=3500/4, length=150)
  ## f0102.dres <- diggle.analysis(f0102.1, f0102.2, plot=T,steps=f.steps)
  
  ## f0103.1 <- matrix(scan("~/mosaics/ferret-th/data/f0103i.txt"),
  ##                   ncol=2, byrow=T)
  ## f0103.2 <- matrix(scan("~/mosaics/ferret-th/data/f0103g.txt"),
  ##                   ncol=2, byrow=T)
  ## f.steps <- seq(from=0, to=3500/4, length=150)
  ## f0103.dres <- diggle.analysis(f0103.1, f0103.2, plot=T,steps=f.steps)

  dev.off()

  ## Save the results of Diggle analysis to disk, for later plotting
  ## with ./plot_fth_diggle.R
  save(f9942.dres, f0102.dres, f0103.dres, file="thdiggle.Rda")
}





## CHAT neurons next?
require(spatstat)
data(amacrine)
scale = 662                             #convert to um.
chat = cbind(amacrine$x, amacrine$y) * scale
chatl = amacrine$marks == "on"
chat.on = chat[chatl,]
chat.of = chat[!chatl,]
plot(chat.on,asp=1)
points(chat.of, pch=19)
chat.w = c(0, 1060, 0, 662)
rect(chat.w[1], chat.w[3], chat.w[2], chat.w[4])


chat.poly = bboxx(bbox(as.points(chat)))
par(mfrow=c(2,3))
chat.steps <- seq(from=0, to=660/4, length=150)
chat.dres <- diggle.analysis(chat.on, chat.of, chat.steps, chat.poly,
                             controls=TRUE,
                             do.csr = FALSE,
                             show.fields = TRUE, plot=T)
Kenv.tor(chat.on, chat.of, chat.poly, s=chat.steps)



######################################################################
## Old code below.
myKenv.label2 <- function(pts1, pts2, poly, nsim, s, quiet=FALSE,
                          plot=FALSE, need.ri=FALSE) {
  ## This function based on Kenv.label from splancs library.
  ## This adapted version also computes the regularity index.

  getri <- function(pts, poly) {
    this.bbox <- bbox(poly)
    vorcr(pts[,1], pts[,2],
          min(this.bbox[,1]), max(this.bbox[,1]),
          min(this.bbox[,2]), max(this.bbox[,2]))$cr
  }
  
  k1.hat <- khat(pts1, poly, s)
  k2.hat <- khat(pts2, poly, s)
  k12.hat<- k12hat(pts1, pts2, poly, s)
  kdiff.real <- k1.hat - k2.hat
  u <- rep(0, length=nsim+1)

  kres <- matrix(0, nrow=nsim+1, ncol=length(kdiff.real))
  kres[1,] <- kdiff.real


  if (need.ri) {
    ri1 <- ri2 <- real(length=nsim+1)
    ri1[1] <- getri(pts1, poly)
    ri2[1] <- getri(pts2, poly)
  } else {
    ri1 <- ri2 <- NULL
  }

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
    
    diffk <- k1sim - k2sim
    kmax <- pmax(kmax, diffk)
    kmin <- pmin(kmin, diffk)

    kres[isim+1,] <- diffk
    if (need.ri) {
      ri1[isim+1] <- getri(labpts[[1]], poly)
      ri2[isim+1] <- getri(labpts[[2]], poly)
    }
  }

  if (plot) {
    plot(s, kdiff.real, type="l", col="red", lty=2,
         xlab=expression(paste("distance (", mu, "m)")), 
         ylab = "L[12]",
         ylim = c( min(c(kmin, kdiff.real)), max(c(kmax, kdiff.real))))
    lines(s, kmin); lines(s, kmax)
  }
  
  list(lower = kmin, upper = kmax, real=kdiff.real, u=u, rank=rank(u)[1],
       kres=kres, ri1=ri1, ri2=ri2)
}
myKenv.label <- function(pts1, pts2, poly, nsim, s, quiet=FALSE, plot=FALSE) {
  ## This function based on Kenv.label from splancs library.

  
  score.var <- function(k1, k2, k12) {
    ## Score taken from (Diggle, 1986, bottom p122)
    sum(apply( rbind( sqrt(k1), sqrt(k2), sqrt(k12)), 2, var))
  }

  k1.hat <- khat(pts1, poly, s)
  k2.hat <- khat(pts2, poly, s)
  k12.hat<- k12hat(pts1, pts2, poly, s)
  kdiff.real <- k1.hat - k2.hat
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
    
    diffk <- k1sim - k2sim
    kmax <- pmax(kmax, diffk)
    kmin <- pmin(kmin, diffk)
  }

  if (plot) {
    plot(s, kdiff.real, type="l", col="red", lty=2,
         xlab=expression(paste("distance (", mu, "m)")), 
         ylab=expression(L[12]),
         ylim = c( min(c(kmin, kdiff.real)), max(c(kmax, kdiff.real))))
    lines(s, kmin); lines(s, kmax)
  }
  
  list(lower = kmin, upper = kmax, real=kdiff.real, u=u, rank=rank(u)[1])
}
