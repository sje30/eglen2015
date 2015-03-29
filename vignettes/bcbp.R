## Plot the BCBP data sets.
## 2015-03-24

require(splancs)

## ---- bcbp-compute
bcbp.dat.file = system.file("extdata/bcbp2/bc_bcbp_f1.dat", package="eglen2015")
bcbp.w.file   = system.file("extdata/bcbp2/bc_bcbp_f1.w",   package="eglen2015")

bcbp = as.matrix(read.table(bcbp.dat.file))
bcbp.w = as.matrix(read.table(bcbp.w.file))

## type 1 = blue cone (BC)
## type 2 = blue cone bipolar (BCBP)

bcbp.computek12 = function(smax) {
  id1 = which(bcbp[,3]==1)
  pts1 = bcbp[ id1, 1:2]
  pts2 = bcbp[-id1, 1:2]
  poly = bboxx(bbox(matrix(bcbp.w,2,2)))
  steps = 0:smax
  k12.hat = k12hat(pts1, pts2, poly, steps)
  cbind(steps, sqrt(k12.hat/pi))
}

## ---- bcbp-draw

pdf(file='bcbp-fig.pdf', width=6.8, height=3.5)
par(mfrow=c(1,2), mar=c(3,1, 0.1, 1.5), mgp=c(2, 0.8, 0))
plot(NA, asp=1, xaxs='i', yaxs='i', bty='n',
     xlim=bcbp.w[1:2], ylim=bcbp.w[3:4],
     xlab='', ylab='', xaxt='n', yaxt='n')
symbols(bcbp[,1], bcbp[,2],
        circles=ifelse(bcbp[,3]==1, 5, 4), # 10um for BC, 8um for BCBP
        bg=ifelse(bcbp[,3]==1, "blue", "red"),
        inches=FALSE, add=TRUE)
rect(bcbp.w[1], bcbp.w[3], bcbp.w[2], bcbp.w[4])
segments(50, 30, 150, lend="butt", lwd=3,xpd=NA)
legend(x=220, y=20, legend=c('BC', 'BCBP'),
       horiz=TRUE, xpd=NA, pch=19, col=c("blue", "red"))
mtext(side=3, at=0, 'A')

## Now do the K function.
par(mar=c(4,3, 0.1, 0.5))
smax = 60
bcbpf1.k12 = bcbp.computek12(smax)
plot(bcbpf1.k12, type='l', asp=1, xlim=c(0,smax), ylim=c(0, smax),
     bty='n', col='red', las=1,
     xlab=expression('distance (' * mu * 'm)'),
     ylab=expression(L[12]))
segments(0, 0, smax, smax, lty=2)
legend(5, 60, legend=c('data', 'CSR'),
       horiz=FALSE,
       lty=c(1,2), col=c('red', 'black'))
mtext(side=3, at=-97, 'A', cex=1.25, line=-1,xpd=NA)
mtext(side=3, at=-15, 'B', cex=1.25, line=-1,xpd=NA)
dev.off()



