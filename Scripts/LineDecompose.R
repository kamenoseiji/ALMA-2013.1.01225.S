#-------- Constants
sigma2FWHM <- 2.0*sqrt(2* log(2))
gaussInteg <- sqrt(2* pi)

#-------- Functions
Gauss2FWHM <- function(sigma, sigmaErr, Tau, TauErr){
	FWHM <-  sigma2FWHM* sigma
	FWHMerr <- sigma2FWHM* sigmaErr
	IntegTau <- Tau* sigma* sqrt(2* pi)
	IntegTauErr <- IntegTau* sqrt((sigmaErr / sigma)^2 + (TauErr / Tau)^2 )
	text_sd <- sprintf('FWHM=%.1f +- %.1f  Peak=%.4f +- %.4f IntegTau=%.1f +- %.1f\n', FWHM, FWHMerr, Tau, TauErr, IntegTau, IntegTauErr)
	cat(text_sd)
}

GaussPlot <- function(a, b, c, lineType=2, color='red'){
	veloc <- seq(1000, 2000, by=1)
	Tau <- a* exp( -0.5*((veloc - b) / c)^2)
	lines(veloc, Tau, lty=lineType, col=color)
	return()
}

setwd('~/Projects/NGC1052_work')

#fileName <- 'CO_J=2-1.data'
#fileName <- 'CO_J=3-2.data' # 2 Gaussian
#fileName <- 'HCN_J=4-3_v=0.data'
#fileName <- '34SO_v=0_6(5)-5(4).data'
#fileName <- 'SO_7_8-7_7.data'
#fileName <- 'SO_5_5-4_4.data'
#fileName <- 'SO_J_N=7_8-6_7.data'
#fileName <- 'SO_J_N=8_8-7_7.data'  # 2 Gaussian
#fileName <- 'SO_J_N=9_8-8_7.data'
#fileName <- 'SO2_v=0_16(3,13)-16(2,14).data'
#fileName <- 'SO2_v=0_22(2,20)-22(1,21).data'
#fileName <- 'SO2_v=0_25(3,23)-25(2,24).data'
#fileName <- 'H2O_5(5,0)-6(4,3)_v2=1.data'
#fileName <- 'CS_J=7-6.data'
#fileName <- 'HCN_J=4-3_v=0.data'
#fileName <- 'HCO+_v2=1_J=4-3_l=1f.data'
#fileName <- 'HCO+_v2=1_J=4-3_l=1e.data'  # 3 Gaussians
#fileName <- 'SO2_v=0_40(4,36)-40(3,37).data' # 3 Gaussians
#fileName <- 'SO2_v=0_9(4,6)-9(3,7).data'  # single
#fileName <- 'CN_N=3-2_J=7_2-5_2.data' # 2 Gaussian
#fileName <- 'HCO+_v2=1_J=4-3_l=1f.data' # 2 Gaussian
#fileName <- 'SO2_v=0_6(4,2)-6(3,3).data # 2 Gaussian
#fileName <- 'SO2_v=0_15(4,12)-15(3,13).data' # 3 Gaussian

SingleGauss <- T
DoubleGauss <- F
TripleGauss <- F
DF <- read.table(fileName, header=T)
F2Vfit <- lm(formula=VLSR ~ freq, DF)	# frequency -> velocity
V2Ffit <- lm(formula=freq ~ VLSR, DF)	# velocity -> frequency
	
newDF <- data.frame(VLSR = seq(min(DF$VLSR), max(DF$VLSR), by=1))
newDF$freq  <- as.numeric(predict(V2Ffit, newDF))
#-------- 1-line Gaussian
if(SingleGauss){
	fit1 <- nls(formula=Tau~Tau0* exp(-0.5*((VLSR - V0)/sigma)^2), data=DF, start=list(Tau0=max(DF$Tau), V0=1492, sigma=70))
	newDF$Tau <- predict(fit1, newDF)
}
DF$resid <- DF$Tau - as.numeric(predict(fit1, DF))

#-------- 2-line Gaussian
if(DoubleGauss){
	#fit2 <- nls(formula=Tau~Tau0* exp(-0.5*((VLSR - V0)/sigma0)^2) + Tau1* exp(-0.5*((VLSR - V1)/sigma1)^2), data=DF, start=list(Tau0=0.14, V0=1500, sigma0=60, Tau1=0.06, V1=1900, sigma1=60)) # for CO(3-2) and H13CN
	# fit2 <- nls(formula=Tau~Tau0* exp(-0.5*((VLSR - V0)/fixedSigma)^2) + Tau1* exp(-0.5*((VLSR - V1)/fixedSigma)^2), data=DF, start=list(Tau0=0.025, V0=1450, Tau1=0.015, V1=1700, fixedSigma=60)) # for CN J=7/2-5/2 and J=5/2-3/2
	# fit2 <- nls(formula=Tau~Tau0* exp(-0.5*((VLSR - V0)/fixedSigma)^2) + Tau1* exp(-0.5*((VLSR - V1)/fixedSigma)^2), data=DF, start=list(Tau0=0.06, V0=1480, Tau1=0.01, V1=1630, fixedSigma=60)) # for SO (8,8)-(7,7) and HC15N
	# fit2 <- nls(formula=Tau~Tau0* exp(-0.5*((VLSR - V0)/fixedSigma)^2) + Tau1* exp(-0.5*((VLSR - V1)/fixedSigma)^2), data=DF, start=list(Tau0=0.032, V0=1460, Tau1=0.03, V1=1740, fixedSigma=70)) # for HCO+ v2=1 l=1f and SO2 17(4,14) - 17(3,15)
	#fit2 <- nls(formula=Tau~Tau0* exp(-0.5*((VLSR - V0)/fixedSigma)^2) + Tau1* exp(-0.5*((VLSR - V1)/fixedSigma)^2), data=DF, start=list(Tau0=0.015, V0=1550, Tau1=0.005, V1=1450, fixedSigma=70)) # for SO2 15(4,12) - 15(3,13), 13(4,10) - 13(13,11), and 11(4,8) - 11(3,9)
	newDF$Tau <- predict(fit2, newDF)
}

#-------- 3-line Gaussian
if(TripleGauss){
	fit3 <- nls(formula=Tau~Tau0* exp(-0.5*((VLSR - V0)/fixedSigma)^2) + Tau1* exp(-0.5*((VLSR - V1)/fixedSigma)^2) + Tau2* exp(-0.5*((VLSR - V2)/fixedSigma)^2), data=DF, start=list(Tau0=0.1, V0=1340, Tau1=0.035, V1=1490, Tau2=0.77, V2=1750, fixedSigma=63)) # for HCO+_v2=1_l=1e, HCO+v=0, and HCN_v2=1_l=1f
	# fit3 <- nls(formula=Tau~Tau0* exp(-0.5*((VLSR - V0)/fixedSigma)^2) + Tau1* exp(-0.5*((VLSR - V1)/fixedSigma)^2) + Tau2* exp(-0.5*((VLSR - V2)/fixedSigma)^2), data=DF, start=list(Tau0=0.01, V0=1270, Tau1=0.001, V1=1460, Tau2=0.09, V2=1560, fixedSigma=60))
	#fit3 <- nls(formula=Tau~Tau0* exp(-0.5*((VLSR - V0)/fixedSigma)^2) + Tau1* exp(-0.5*((VLSR - V1)/fixedSigma)^2) + Tau2* exp(-0.5*((VLSR - V2)/fixedSigma)^2), data=DF, start=list(Tau0=0.012, V0=1560, Tau1=0.004, V1=1470, Tau2=0.003, V2=1390, fixedSigma=60)) # for SO2 15(4,12) - 15(3,13), 13(4,10) - 13(13,11), and 11(4,8) - 11(3,9)
	newDF$Tau <- predict(fit3, newDF)
}

plot(newDF$VLSR, newDF$Tau, xlab='VLSR (km/s)', ylab='Optical depth', main=fileName, type='l')
points(DF$VLSR, DF$Tau, pch=20)

#---- Single-Gaussian
if(SingleGauss){
	cat(sprintf('VLSR = %.1f +- %.1f weightedVLSR = %.1f ', as.numeric(coef(fit1)['V0']), as.numeric(summary(fit1)$parameters[,2]['V0']), sum(DF$VLSR* DF$Tau)/sum(DF$Tau)))
	Gauss2FWHM(as.numeric(coef(fit1)['sigma']), as.numeric(summary(fit1)$parameters[,2]['sigma']), as.numeric(coef(fit1)['Tau0']), as.numeric(summary(fit1)$parameters[,2]['Tau0']))
}

#-------- For decomposition at 356 GHz
if(TripleGauss){
	fixedSigma <- coef(fit3)['fixedSigma']
	GaussPlot(coef(fit3)['Tau0'], coef(fit3)['V0'], fixedSigma, 2, 'blue')
	GaussPlot(coef(fit3)['Tau1'], coef(fit3)['V1'], fixedSigma, 2, 'green')
	GaussPlot(coef(fit3)['Tau2'], coef(fit3)['V2'], fixedSigma, 2, 'orange')

	F2Vfit <- lm(formula=VLSR ~ freq, DF)
	V2Ffit <- lm(formula=freq ~ VLSR, DF)

	freqRef <- 356.549   # HCO+_v2=1_J=4-3_l=1e
	#freqRef <- 341.403   # SO2 v=0 40(4,36)-40(3,37)
	#---- 
	# freqRest<- 341.67395780 # SO2 v=0 36(5,31)-36(4,32)', restfreq=341.67395780
	freqRest<- 356.73424 # HCO+ J=4-3', restfreq=356.73424
	FreqTopo <- as.numeric(predict(V2Ffit, data.frame(VLSR=coef(fit3)['V0'])))
	LineVeloc <- as.numeric(predict(F2Vfit, data.frame(freq=FreqTopo - freqRest + freqRef)))
	cat(sprintf('VLSR = %.1f +- %.1f ', LineVeloc, as.numeric(summary(fit3)$parameters[,2]['V0'])))
	Gauss2FWHM(fixedSigma, as.numeric(summary(fit3)$parameters[,2]['fixedSigma']), as.numeric(coef(fit3)['Tau0']), as.numeric(summary(fit3)$parameters[,2]['Tau0']))

	#---- 
	# freqRest<- 341.403 # 'SO2 v=0 40(4,36)-40(3,37)', restfreq=341.403
	freqRest<- 356.549 # HCO+ v2=1 J=4-3 l=1e', restfreq=356.549
	FreqTopo <- as.numeric(predict(V2Ffit, data.frame(VLSR=coef(fit3)['V1'])))
	LineVeloc <- as.numeric(predict(F2Vfit, data.frame(freq=FreqTopo - freqRest + freqRef)))
	cat(sprintf('VLSR = %.1f +- %.1f ', LineVeloc, as.numeric(summary(fit3)$parameters[,2]['V1'])))
	Gauss2FWHM(fixedSigma, as.numeric(summary(fit3)$parameters[,2]['fixedSigma']), as.numeric(coef(fit3)['Tau1']), as.numeric(summary(fit3)$parameters[,2]['Tau1']))

	#---- 
	# freqRest<- 341.27552080 # SO2 v=0 21(8,14)-22(7,15)', restfreq=341.27552080
	freqRest<- 356.256 # HCN v2=1 J=4-3 l=1f', restfreq=356.256
	FreqTopo <- as.numeric(predict(V2Ffit, data.frame(VLSR=coef(fit3)['V2'])))
	LineVeloc <- as.numeric(predict(F2Vfit, data.frame(freq=FreqTopo - freqRest + freqRef)))
	cat(sprintf('VLSR = %.1f +- %.1f ', LineVeloc, as.numeric(summary(fit3)$parameters[,2]['V2'])))
	Gauss2FWHM(fixedSigma, as.numeric(summary(fit3)$parameters[,2]['fixedSigma']), as.numeric(coef(fit3)['Tau2']), as.numeric(summary(fit3)$parameters[,2]['Tau2']))
}

#-------- For CO(3-2) + H13CN
if(DoubleGauss){
	fixedSigma <- coef(fit2)['fixedSigma']
	GaussPlot(coef(fit2)['Tau0'], coef(fit2)['V0'], fixedSigma, 2, 'blue')
	GaussPlot(coef(fit2)['Tau1'], coef(fit2)['V1'], fixedSigma, 2, 'red')

	F2Vfit <- lm(formula=VLSR ~ freq, DF)
	V2Ffit <- lm(formula=freq ~ VLSR, DF)

	#freqRef <- 340.247874   # CN N=3-2 J=7/2-5/2
	#freqRef <- 344.3106     # SO J_N=8_8-7_7', restfreq=344.3106
	#freqRef <- 358.242     # HCO+ v2=1 J=4-3 l=1f', restfreq=358.242
	freqRef <- 357.241     # SO2 v=0 15(4,12)-15(3,13)', restfreq=357.241
	#---- 
	# freqRest<- 340.247874 #CN N=3-2 J=7/2-5/2', restfreq=340.247874
	#freqRest<- 344.3106 # SO J_N=8_8-7_7', restfreq=344.3106
	# freqRest<- 358.242 # HCO+ v2=1 J=4-3 l=1f', restfreq=358.242
	freqRest<- 357.165 # 13(4,10)-13(3,11)', restfreq=357.165
	FreqTopo <- as.numeric(predict(V2Ffit, data.frame(VLSR=coef(fit2)['V0'])))
	LineVeloc <- as.numeric(predict(F2Vfit, data.frame(freq=FreqTopo - freqRest + freqRef)))
	cat(sprintf('VLSR = %.1f +- %.1f ', LineVeloc, as.numeric(summary(fit2)$parameters[,2]['V0'])))
	Gauss2FWHM(coef(fit2)['fixedSigma'], as.numeric(summary(fit2)$parameters[,2]['fixedSigma']), as.numeric(coef(fit2)['Tau0']), as.numeric(summary(fit2)$parameters[,2]['Tau0']))

	#---- 
	# freqRest<- 340.031567 # CN N=3-2 J=5/2-3/2', restfreq=340.031567
	# freqRest<- 344.20031990 # HC15N v=0 J=4-3', restfreq=344.20031990
	# freqRest<- 357.963 # SO2 v=0 17(4,14)-17(3,15)', restfreq=357.963
	freqRest<- 357.388 # 11(4,8)-11(3,9)', restfreq=357.388
	FreqTopo <- as.numeric(predict(V2Ffit, data.frame(VLSR=coef(fit2)['V1'])))
	LineVeloc <- as.numeric(predict(F2Vfit, data.frame(freq=FreqTopo - freqRest + freqRef)))
	cat(sprintf('VLSR = %.1f +- %.1f ', LineVeloc, as.numeric(summary(fit2)$parameters[,2]['V1'])))
	Gauss2FWHM(coef(fit2)['fixedSigma'], as.numeric(summary(fit2)$parameters[,2]['fixedSigma']), as.numeric(coef(fit2)['Tau1']), as.numeric(summary(fit2)$parameters[,2]['Tau1']))
}