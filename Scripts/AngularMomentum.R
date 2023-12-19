AM <- read.csv('~/Projects/NGC1052_work/AngularMomentum.csv')
AM$pc <- AM$Radius * 85	# conversion arcsec -> pc
AM$am <- abs(AM$VLSR - 1492.0)* AM$pc # angular momentum

crossRadius <- (AM[AM$MC == 'Torus',]$am/sqrt(AM[AM$MC == 'Torus',]$pc) / AM[AM$MC == 'CND',]$am * AM[AM$MC == 'CND',]$pc)^2	# 14.395 pc
crossFact1 <- AM[AM$MC == 'CND',]$am / AM[AM$MC == 'CND',]$pc		# 245 km/s
crossFact2 <- AM[AM$MC == 'Torus',]$am / sqrt(AM[AM$MC == 'Torus',]$pc)	# 929.516 km/s/pc^1/2

pdf('AngularMomentumNGC1052.pdf')
plot(AM$pc, AM$am, log='xy', pch=20, xlab='Distance from the nucleus [pc]', ylab='Specific angular momentum [pc km/s]', main='Molecular gas in NGC1052', xlim=c(2, 1000), asp=1); text(AM$pc, AM$am, labels=AM$MC, adj=c(-0.1,1), cex=0.7)
lines(c(crossRadius, 500), c(crossRadius*245, 500*crossFact1))
lines(c(5,crossRadius), c(5*crossFact1, crossRadius*crossFact1), lty=2)
lines(c(2, crossRadius), c(crossFact2*sqrt(2), crossFact2*sqrt(crossRadius)))
lines(c(crossRadius, 100), c(crossFact2*sqrt(crossRadius),crossFact2*sqrt(100)), lty=2)
dev.off()