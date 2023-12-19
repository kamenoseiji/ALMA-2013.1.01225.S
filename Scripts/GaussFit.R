library(xtable)

fitLine <- function(LDF){	# Label, param, offset=0, velRange = c(1200, 1800)){
	File <- sprintf('%s.data', LDF$label)
	DF <- read.table(File, header=T)
	DF <- DF[((DF$VLSR > LDF$VB) & (DF$VLSR < LDF$VE)),]
	DF$Tau <- DF$Tau - LDF$offset
	fit <- nls(DF, formula=Tau ~ a1* exp(-0.5*((VLSR - v1)/d1)^2)
							   + a2* exp(-0.5*((VLSR - v2)/d2)^2), start = LDF[c('a1', 'v1', 'd1', 'a2', 'v2', 'd2')],
							   control=nls.control(maxiter=1000))
    pdf(sprintf('%s.fit.pdf', LDF$label))
    plotVeloc <- seq(1200, 1800, by=10)
	plot(DF$VLSR, DF$Tau, pch=20, cex=0.5, type='n', xlab='LSR Velocity [km/s]', ylab='Optical Depth', main=LDF$label )
	lines(plotVeloc, predict(fit, list(VLSR=plotVeloc)), col='red')
	points(DF$VLSR, DF$Tau, pch=20, cex=0.5)
	dev.off()
	return( as.numeric(coef(fit)) )
}

integLine <- function(params){
	amp <- params[seq(1, length(params), by=3)]
	vel <- params[seq(2, length(params), by=3)]
	wid <- abs(params[seq(3, length(params), by=3)])
	return( sqrt(2.0*pi) * amp %*% wid )
}

lineLists <- c('CO_J=2-1', 'SO_5_5-4_4', 'SO_7_8-7_7', 'HCO+_J=4-3', 'HCN_J=4-3_v=0', 'HCN_v2=1_J=4-3_l=1f', 'CO_J=3-2', 'H13CN_v=0_J=4-3', 'SO_J_N=8_8-7_7', 'SO_J_N=9_8-8_7', 'CS_J=7-6', 'SO_J_N=7_8-6_7')
offsets   <- c(0.0, 0.0, -0.001, 0.0, 0.0, 0.0, 0.005, 0.0, 0.0, 0.002, 0.0, 0.0)
velBegin  <- c(1200, 1300, 1200, 1250, 1150, 1350, 1250, 1380, 1270, 1200, 1200, 1230)
velEnd    <- c(1800, 1800, 1800, 1700, 1800, 1800, 1700, 1800, 1800, 1700, 1770, 1680)
a1 <- c(0.1, 0.03, 0.006, 0.1, 0.15, 0.08, 0.14, 0.06, 0.06, 0.056, 0.03, 0.045)
a2 <- c(0.01, 0.001, 0.001, 0.04, 0.02, 0.01, 0.01, 0.01, 0.01, 0.004, 0.005, 0.01)
v1 <- c(1500, 1491, 1500, 1480, 1480, 1490, 1490, 1490, 1490, 1484, 1500, 1470)
v2 <- c(1650, 1680, 1350, 1630, 1550, 1670, 1530, 1660, 1380, 1350, 1300, 1590)
d1 <- c(70, 63, 30, 50, 60, 60, 50, 70, 70, 71, 80, 50)
d2 <- c(100, 57, 90, 80, 150, 40, 100, 20, 40, 29, 80, 20)
lineDF <- data.frame(label=lineLists, offset=offsets, VB=velBegin, VE=velEnd, a1, a2, v1, v2, d1, d2)
integTau <- numeric(0)
#---- Each line
for(index in 1:nrow(lineDF)){
	params <- fitLine(lineDF[index,])
	lineDF[index,][c('a1', 'v1', 'd1', 'a2', 'v2', 'd2')] <- params
	integTau[index] <- integLine(params)
}
lineDF$d1 <- abs(lineDF$d1)
lineDF$d2 <- abs(lineDF$d2)
lineDF$EW <- integTau
print(lineDF)