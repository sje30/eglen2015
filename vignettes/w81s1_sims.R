## ---- get-rgc-data
##rgc.of = read.table('~/mosaics/data/w81s1of.txt')
##rgc.on = read.table('~/mosaics/data/w81s1on.txt')
rgc.of.file = system.file("extdata/w81s1/w81s1of.txt", package="eglen2015")
rgc.on.file = system.file("extdata/w81s1/w81s1on.txt", package="eglen2015")
rgc.w.file  = system.file("extdata/w81s1/w81s1w.txt",  package="eglen2015")

rgc.of = read.table(rgc.of.file)
rgc.on = read.table(rgc.on.file)
rgc.w  = scan(rgc.w.file)

## ---- show-univ
rgc.soma.rad = 8
par(bty='n', mfrow=c(1,2))
plot(rgc.on, asp=1, type='n', xlab='', ylab='')
symbols(rgc.on[,1], rgc.on[,2],
asp=1, circles=rep(rgc.soma.rad, nrow(rgc.on)),
        inches=FALSE, add=TRUE)
rect(rgc.w[1], rgc.w[3], rgc.w[2], rgc.w[4])



