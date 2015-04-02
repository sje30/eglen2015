## chickcones.R --- show the Kram et al (2010) chick photoreceptors.
## ---- plot-chickcones
kram.dn4.file = system.file("extdata/kram/dn4.csv", package="eglen2015")
kram.data = read.csv(kram.dn4.file, skip=1)
cols = scan(kram.dn4.file, what="character", nlines=1, sep=',')
cols[9] = "black"
xcol = 1;

pdf(file="chickcones.pdf", width=4, height=4)
par(mar=rep(0.5, 4))
plot(NA, xlim=c(100,200), ylim=c(100,200), asp=1, xaxt='n', yaxt='n')
for (i in 1:5) {
  x = kram.data[,xcol]; y=kram.data[,xcol+1]
  col = cols[xcol]
  ##points(x, y, pch=19, col=col, cex=0.3)
  symbols(x, y, circles=rep(1.1, length(x)), lwd=0.01,
          inches=FALSE, add=TRUE, bg=col, cex=0.3)
  xcol = xcol + 2                       #get ready for next type
}
segments(100, 95, 120, lend="butt", lwd=4,xpd=NA)
dev.off()
