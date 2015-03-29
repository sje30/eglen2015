## Adapted from: ~/mosaics/beta_rgc/dmin/dmin2a.R
## 2015-02-07
## Analyse and model beta maps.
## This script should run by itself.

## If a chunk does not have a name, then if there is an empty chunkname
## read in, it will be included! (Bad explanation).  Giving the preamble above
## a chunkname gets around this.

## ---- get-rgc-data
## load in data files.
##rgc.of = read.table('~/mosaics/data/w81s1of.txt')
##rgc.on = read.table('~/mosaics/data/w81s1on.txt')
rgc.of.file = system.file("extdata/w81s1/w81s1of.txt", package="eglen2015")
rgc.on.file = system.file("extdata/w81s1/w81s1on.txt", package="eglen2015")
rgc.w.file  = system.file("extdata/w81s1/w81s1w.txt",  package="eglen2015")
rgc.of = as.matrix(read.table(rgc.of.file))
rgc.on = as.matrix(read.table(rgc.on.file))
rgc.w  = scan(rgc.w.file)


## ---- setup-betamap

require(sjedmin)
require(sjedist)
library(spatstat)

## TODO attach("~/mosaics/data/bivariate_mosaics.Rda")


## First set up h for all sims.
hpar <- function(d,theta) {
  ## Choice of h() suggested by Peter.
  delta<-theta[1]
  phi<-theta[2]
  alpha<-theta[3]
  res <- (0*(d<delta)) + (d>=delta)*(1-exp(-((d-delta)/phi)^alpha))
  if (any (is.nan(res)))
    res[ which(is.nan(res))] <- 0

  res
}

ranking <- function(arr) {
  ## Evaluate fit of row 1 (the real data) with remaining rows (simulations)
  ## using equations 6.3--6.5 from (Diggle, 2002, page 89).
  n.s <- nrow(arr)
  u <- rep(0, n.s)
  for (i in 1:n.s) {
    ave.i <- apply( arr[-i,], 2, sum) / (n.s - 1)
    u[i] <- sum((ave.i - arr[i,])^2)
    ##u[i] <- max( abs( ave.i - arr[i,]) )
  }
  signif( (rank(-u)[1])/n.s, 5)
}

ranking.k <- function(arr) {
  ## Assume we have a K function; cf it against the pi t^2 version.
  n.s <- nrow(arr)
  u <- rep(0, n.s)
  for (i in 1:n.s) {
    ave.i <- apply( arr[-i,], 2, sum) / (n.s - 1)
    u[i] <- sum((ave.i - arr[i,])^2)
    ##u[i] <- max( abs( ave.i - arr[i,]) )
  }
  signif( (rank(-u)[1])/n.s, 5)
}


make.ci <- function(da, fns=c("l1", "f1")) {
  ## Make the confidence intervals for each of the arrays mentioned in fns.
  ## Return these confidence intervals as a list.

  res <- list()
  for (f in fns) {
    arr <- da$get()[[f]]
    p <- 0.05                             # 5% for 95% confidence intervals 
    ci <- apply(arr, 2, quantile, probs=c(p/2, 1-(p/2)))
    
    ## This would be nice, but doesn't work.
    ## attr(dist.arr$get()[[f]], "ci") <- ci
    
    new <- list(ci); names(new) <- f
    res <- c(res, new)
    
    res
  }
}


kranking <- function(arr) {
  ## Assume we have a L function; cf it against the pi t^2 version.
  n.s <- nrow(arr)
  u <- rep(0, n.s)

  t <- as.real(colnames(arr))
  theor <- sqrt(pi *t^2)
  for (i in 1:n.s) {
    diffs <- ( (arr[i,]) - theor)^2
    u[i] <- sum( diffs)
  }
  print(u)
  signif( (rank(-u)[1])/n.s, 5)
}

lplot <- function(arr, ci=NULL, r=NULL, ylab, do.krank=FALSE, ...) {
  ## Plot a spatial distribution (K, F, G).
  ## Real data is shown in red; simulation envelope in black.
  if (is.null(r)) {
    r <- colnames(arr)
  }
  if (missing (ylab)) {
    ylab <- attributes(arr)$name
    if (is.null(ylab))
      ylab <- deparse(substitute(arr))
  }
  ## are we plotting a g function
  gfunction <- (max(arr[1,]) < 1.1)
  if (gfunction) {
    ylim <- c(0,1)
  } else {
    ylim <- c(0, 150)
  }

  if (do.krank) {
    pval <- kranking(arr)
  } else {
    pval <- ranking(arr)
  }
  plot(r, arr[1,], col='black', type='l', bty='n',
       ##main=ranking(arr),
       ylim=ylim,
       ylab=ylab, xlab='')
  title(xlab=expression(paste("distance (", mu, "m)")), mgp=c(2,1,0))
  ##text(0, 140, paste('p = ', pval), adj=0)
  mtext(paste('p = ', pval), adj=0.1, side=3, line=-1, cex=0.8)
  ##title(main= paste("p =", ranking(arr)), line=-1)
  if (is.null(ci)) {
    ## if no CI given, just plot envelope.
    min.line <- apply(arr[-1,], 2, min)
    max.line <- apply(arr[-1,], 2, max)
    lines(r, min.line, lty=1);   lines(r, max.line, lty=1)
  } else {
    ## We have the CI data, so plot it.
    lines(r, ci[1,], lty=2)
    lines(r, ci[2,], lty=2)
  }
}

w81s1.plot <- function(){
  ps.options(family="Helvetica")
  postscript(file="|psfbb>w81s1_gof2t.ps",
             horiz=F, onefile=FALSE, width=7, height=9)
  par(mar=c(3,4, 1, 1), mfrow=c(4,2), bty='n', las=1,oma=c(1,0,0,0))
  lplot(w81s1.arr$get()$l1,  w81s1.ci$l1,
        ylab=expression(L[1]))
  lplot(w81s1.arr$get()$l2,  w81s1.ci$l2,
        ylab=expression(L[2]))
  lplot(w81s1.arr$get()$l0,  w81s1.ci$l0,
        ylab=expression(L[1+2]))
  lplot(w81s1.arr$get()$l12, w81s1.ci$l12,
        ylab=expression(L[12]))
  lplot(w81s1.arr$get()$g1,  w81s1.ci$g1,  ylab=expression(G[1]))
  lplot(w81s1.arr$get()$g2,  w81s1.ci$g2,  ylab=expression(G[2]))
  plot.spat.opp2(w81s1.arr$get()$opp, cex=0.1, real.col='black')
  plot.spat.ri3.v2(w81s1.arr$get()$ri3, cex=0.1, real.col='black')
  ##mtext('W81s1 29 Oct 2004; bivariate pipp theta1=theta2=(45, 41.53, 2.66), d12=18', side=1, outer=T)
  dev.off()
}

m623.plot <- function() {
  ps.options(family="Helvetica")
  postscript(file="|psfbb>m623_gof2t.ps", horiz=F,
             onefile=FALSE, width=7, height=9)
  par(mar=c(3,4, 1, 1), mfrow=c(4,2), bty='n', las=1,oma=c(1,0,0,0))
  lplot(m623.arr$get()$l1,  m623.ci$l1,   ylab=expression(L[1]))
  lplot(m623.arr$get()$l2,  m623.ci$l2,   ylab=expression(L[2]))
  lplot(m623.arr$get()$l0,  m623.ci$l0,   ylab=expression(L[1+2]))
  lplot(m623.arr$get()$l12, m623.ci$l12, ylab=expression(L[12]))
  lplot(m623.arr$get()$g1,   m623.ci$g1,  ylab=expression(G[1]))
  lplot(m623.arr$get()$g2,   m623.ci$g2,  ylab=expression(G[2]))
  plot.spat.opp2(m623.arr$get()$opp, cex=.1, real.col='black')
  plot.spat.ri3.v2(m623.arr$get()$ri3, cex=.1, real.col='black')
  ##mtext('W81s1 29 Oct 2004; bivariate pipp theta1=theta2=(45, 41.53, 2.66), d12=18', side=1, outer=T)
  dev.off()
}

m623.handplot <- function() {
  ps.options(family="Helvetica")
  postscript(file="|psfbb>m623_handgof.ps", horiz=F,
             onefile=FALSE, width=7, height=7)
  par(mar=c(3,4, 1, 1), mfrow=c(2,2), bty='n', las=1,oma=c(1,0,0,0))
  lplot(m623.arr$get()$l1,  m623.ci$l1,   ylab=expression(L[1]))
  lplot(m623.arr$get()$l12, m623.ci$l12, ylab=expression(L[12]))
  lplot(m623.arr$get()$g1,   m623.ci$g1,  ylab=expression(G[1]))
  plot.spat.ri3.v2(m623.arr$get()$ri3, cex=.1, real.col='black')
  dev.off()
}

plot.spat.opp2 <- function(arr,cex=0.5, real.col='black', ...) {
  ## Plot the fraction of opposites

  opp <- arr[-1,]
  boxplot( list(opp[,1], opp[,2], opp[,3], opp[,4]), xaxt ='n',
          ylab="fraction opposite type",
          pch=19, cex=0.2, pars=list(medlwd=1))

##   stripchart(list(arr[-1,1], arr[-1,2], arr[-1,3], arr[-1,4]),
##              method="jitter", pch=19, vertical=TRUE,
##              ylim=c(min(arr), 1), 
##              ##group.names=c(expression(1^{st}), "2", "3", "all"),
##              group.names=rep("", 4),
##              ylab="fraction opposite type", cex=cex, ...)

  axis(1, at=1:4, labels=c(expression(1^{st}),
                    expression(2^{nd}),
                    expression(3^{rd}),
                    "all"))
  dx <- 0.6; i <- 1:4
  ##segments(i-dx, arr[1,], i+dx, arr[1,], lwd=0.6, col=real.col)
  points(i, arr[1,], pch=18, cex=2)

  ##median.sim <- apply(arr[-1,], 2, median)
  ##segments(i-dx, median.sim, i+dx, median.sim, lwd=0.6, lty=2)

}


plot.spat.ri3.v2 <- function(ri3, cex=0.5, ylim=range(ri3),
                          real.col = 'red', ...) {
  ## Plot the regularity indexes.
  res <- list(on=ri3[-1,1], of=ri3[-1,2],on.off=ri3[-1,3])
  boxplot(res, ylab="regularity index",
          names = c("ON", "OFF", "ON+OFF"),
          pch=19, cex=0.2, pars=list(medlwd=1))
##   stripchart(res, vert=T, pch=19, method="jitter",
##              cex=cex, ylim=ylim,
##              group.names=c("ON", "OFF", "ON+OFF"),
##              main="", ylab="regularity index")
  

  i <- 1:3; dx <- 0.3;
  ##segments(i-dx, ri3[1,], i+dx, ri3[1,], lwd=0.6, col=real.col)
  points(i, ri3[1,], pch=18, cex=2)

  ## median.sim <- apply(ri3[-1,], 2, median)
  ## segments(i-dx, median.sim, i+dx, median.sim, lwd=0.6, lty=2)
  ##legend(x=1, y=3.5, lty=c(1,2),
  ##       legend=c("experimental RI", "median RI of sims"))
  
}

h.nopar <- function(pts, win, rs, plot=FALSE) {
  ## Adapted from ~/mosaics/beta_rgc/pipp/plot_nonparh.R
  ##
  ## Compute the non-parametric estimate of h() for a dataset.
  ## Input:
  ## PTS - (N,2) matrix of data points.
  ## WIN - 4 element window.
  ## NTILES - 2-vector giving number of tiles in (x,y) dimension for
  ## quadrature.
  ## RS - vector of breakpoints where h() is estimated.
  ##
  ## Output: list(x,y) estimating the h() function.
    
  ppp <- as.ppp(pts, as.owin(win))
  qs <- quadscheme(ppp)
  x <- ppm(qs, ~1, PairPiece(r = rs), correction="isotropic")
  h <- summary(x)$interaction$printable

  if (any(is.na(h)))
    h[which(is.na(h))] <- 0
  
  ## return the list of (x,y) points.
  
  ## do we want to take the mid points of RS?  RS are the breaks.
  x <- apply(rbind(rs, c(0, rs[1:(length(rs)-1)])), 2, mean)
  list(x=x, y=h)
}

plot.hnopar <- function(arr, xs, p=0.05) {
  ## Plot the non-parametric estimate of h and confidence intervals.

  ci <- apply(arr[-1,], 2, quantile, probs=c(p/2, 1-(p/2)))
  plot(xs, arr[1,], type='l', col='red', main='',
       xlab=expression(paste("distance (", mu, "m)")),
       ylab="h",
       ylim=c(0,1.8))
  title(deparse(substitute(arr)))
  lines(xs, ci[1,]); lines(xs, ci[2,])
}

######################################################################
## End of functions
######################################################################


######################################################################
## Try W81s GOF
######################################################################

## taken from attach("~/mosaics/data/bivariate_mosaics.Rda")


## ---- w81s1-define-h
w81s1.on.par <- c(15.000000, 67.944693,  7.814627)
w81s1.of.par <- c(15.000000, 66.270221,  5.395813)

h11.x <- seq(from=0, to=160, by=5); h22.x <- h11.x
h11.y <- hpar(h11.x, w81s1.on.par)
h22.y <- hpar(h22.x, w81s1.of.par)

h12.x <- seq(from=17.5, to=18.5, by=.05)
h12.y <- ifelse(h12.x>18.0, 1, 0)

vd.num <- 30                            #given less than 100 real pts
w81s.bdpar <- list( steps=seq(from=1, to=150, length=100),
                   vd0.breaks=seq(from=0, to=12000, len=vd.num),
                   vd1.breaks=seq(from=0, to=20000, len=vd.num),
                   ds0.breaks=seq(from=0, to=200, len=vd.num),
                   ds1.breaks=seq(from=0, to=250, len=vd.num),
                   distribs=list(g0=1, g1=1,g2=1,
                     f0=1, f1=1,f2=1,
                     l0=1, l1=1, l2=1, l12=1,
                     vd0=1, vd1=1, vd2=1,
                     ds0=1, ds1=1, ds2=1,
                     opp=1, ri3=1))

## ---- w81s1-compute-mosaics

n1 <- nrow(rgc.on); n2 <- nrow(rgc.of)

nreps = 9
w81s1.arr <- new.dist.arr( sjespatdists.biv(rgc.on, rgc.of, rgc.w, "note",
                                           param=w81s.bdpar), nreps)

w81s1.rs <- seq(from=10, to=200, by=10)



for (i in 1:nreps) {
  if ((i %%10) == 0)
    cat(paste("iteration", i,"\n"))
  nsweeps <- ifelse(i==1, 10, 10)
  sim.on <- rgc.on; sim.of <- rgc.of
  d <- pipp2.lookup(w=rgc.w,
                    pts1=sim.on, pts2=sim.of,
                    n1=n1, n2=n2,
                    h1=h11.y, d1=h11.x,
                    h2=h22.y, d2=h22.x,
                    h12=h12.y, d12=h12.x, tor=FALSE,
                    nsweeps=nsweeps, verbose=FALSE)

  simpts <- cbind(d$x, d$y)
  sim.on <- simpts[1:(d$n1),]; sim.of <- simpts[-(1:(d$n1)),]
  s <- sjespatdists.biv(sim.on, sim.of, rgc.w, "note", param=w81s.bdpar)
  w81s1.arr$set.row(s, i+1)
}

w81s1.sim.on <- sim.on; w81s1.sim.of <- sim.of

w81s1.ci <- make.ci(w81s1.arr, c("l1", "l2", "l0", "l12", "g1", "g2") )

## Following lines made plots for journal article and saved data.
## w81s1.plot()
## save.image("bpipp_all_date.Rda", compress=TRUE)

## last good data sets: bpipp_w81s1_may20_1.Rda

## bpipp_all_jun2_henv.Rda.gz
## This has output from the 99 simulations of h estimates, useful for
## a figure maybe?

## ---- plot-w81s1-real-sim
rgc.soma.rad = 8
par(mfrow=c(1,2), bty='n')
plot(rgc.on, asp=1, type='n', xlab='', ylab='', xaxt='n', yaxt='n')
symbols(rgc.on[,1], rgc.on[,2], circles=rep(rgc.soma.rad, nrow(rgc.on)),
        inches=FALSE, add=TRUE)
symbols(rgc.of[,1], rgc.of[,2], circles=rep(rgc.soma.rad, nrow(rgc.of)),
        inches=FALSE, add=TRUE,bg='black')
rect(rgc.w[1], rgc.w[3], rgc.w[2], rgc.w[4])
legend(920, 600, legend=c('on-centre', 'off-centre'),
       xpd=NA, pch=c(1, 19),cex=0.7)
mtext(side=3, at=0, 'A')
##
## Now show the simulation ##
plot(rgc.on, asp=1, type='n', xlab='', ylab='', xaxt='n', yaxt='n')
symbols(sim.on[,1], sim.on[,2], circles=rep(rgc.soma.rad, nrow(sim.on)),
        inches=FALSE, add=TRUE)
symbols(sim.of[,1], sim.of[,2], circles=rep(rgc.soma.rad, nrow(sim.of)),
        inches=FALSE, add=TRUE,bg='black')
rect(rgc.w[1], rgc.w[3], rgc.w[2], rgc.w[4])
mtext(side=3, at=0, 'B')
## text(grconvertX(c(0.1, 0.5), from='ndc'), grconvertY(rep(0.9,2), from='ndc'),
##      toupper(letters[1:2]), xpd=NA)




## ---- compute-univariate-betamap
pts = rgc.on
n1 = nrow(pts)
univ.sim = pipp.lookup(rgc.w, pts, n1=n1, h = h11.y, d=h11.x)
            

## ---- plot-univariate-betamap
pdf(file='beta_univ.pdf', width=6.5, height=4)
par(mar=c(0.6,0.01,0.5,0.1))
rgc.soma.rad = 8
par(mfrow=c(1,2), bty='n')
plot(NA, asp=1, xaxs='i', yaxs='i',
     xlim=rgc.w[1:2], ylim=rgc.w[3:4],
     xlab='', ylab='', xaxt='n', yaxt='n')
symbols(rgc.on[,1], rgc.on[,2], circles=rep(rgc.soma.rad, nrow(rgc.on)),
        inches=FALSE, add=TRUE, bg='black')
## symbols(rgc.of[,1], rgc.of[,2], circles=rep(rgc.soma.rad, nrow(rgc.of)),
##         inches=FALSE, add=TRUE,bg='black')
rect(rgc.w[1], rgc.w[3], rgc.w[2], rgc.w[4]) 
mtext(side=3, adj=0, 'A',cex=1.5, line=-1.25)
segments(50, 0, 150, lend="butt", lwd=5,xpd=NA)
##
## Now show the simulation ##
plot(NA, asp=1, xaxs='i', yaxs='i',
     xlim=rgc.w[1:2], ylim=rgc.w[3:4],
     xlab='', ylab='', xaxt='n', yaxt='n')
symbols(univ.sim$x, univ.sim$y, circles=rep(rgc.soma.rad, n1),
        inches=FALSE, add=TRUE, bg='black')
rect(rgc.w[1], rgc.w[3], rgc.w[2], rgc.w[4])
mtext(side=3, adj=0, cex=1.5, 'B', line=-1.25)
dev.off()


## ---- khat-univariate-betamap
require(splancs)
window.to.poly = function(window) {
  bboxx(bbox(matrix(window,2,2)))  
}


lhat = function(pts, steps) {
  poly = window.to.poly(rgc.w)
  k = khat(pts, poly, steps)
  cbind(steps, sqrt(k/pi))
}

steps = seq(from=0, to=200)

pdf(file='khat-univbeta.pdf', width=4, height=4)
par(mgp=c(2,0.8,0), mar=c(3.25,3.25,0.5, 0.5))
plot(lhat(rgc.on, steps), type='l',las=1, bty='n', asp=1, col='red',xlab='')
title(xlab=expression(paste("distance (", mu, "m)")),
      ylab='L')
lines(lhat(cbind(univ.sim$x, univ.sim$y), steps), lty=3, col='blue')
segments(0, 0, max(steps), max(steps), lty=2)
legend(x=0, y=200, ##'topleft',
       legend=c('data', 'model', 'CSR'),
       col=c('red', 'blue', 'black'),
       lty=c(1, 3, 2))
dev.off()
