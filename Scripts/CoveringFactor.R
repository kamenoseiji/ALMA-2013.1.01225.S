cvel <- 299792458
kb   <- 1.38e3  # for Jy unit
distance <- 17.5e6				# angular distance 

obsFluxEm <- function(Fcov, Tau, Scont = 0.444981, Tex, freq=356.734e9, r=7.4, H=1.0){
	velDispersion <- 1 + 1.0 / H^2	# velocity dispersion scaled by absorption line width
	Radius <- r / velDispersion		# Radius in pc
	Omega <- 8.0* (Radius / distance)^2 * H
	return(   Scont* (1.0 - Fcov)	# unabsorbed continuum
			+ Scont* Fcov* exp(-Tau) 	# absorbed continuum
			+ 2.0* kb* Tex* (1.0 - exp(-Tau / velDispersion ))* (freq / cvel)^2 * Fcov* Omega )	# emission
}

FluxEm <- function(Fcov, Tau, Tex, freq=356.734e9, r=7.4, H=1.0){
	velDispersion <- 1 + 1.0 / H^2	# velocity dispersion scaled by absorption line width
	Radius <- r / velDispersion		# Radius in pc
	Omega <- 8.0* (Radius / distance)^2 * H
	return(2.0* kb* Tex* (1.0 - exp(-Tau / velDispersion ))* (freq / cvel)^2 * Fcov* Omega )
}


fcov <- seq(0.1, 1.0, by=0.001)
tauH12CN <- seq(2.0, 60, by=0.1)
tauH13CN <- tauH12CN / 50
fluxSigmaH12CN <- 0.0003277
fluxSigmaH13CN <- 0.0002772
fluxSigmaCO21  <- 0.0002739
#HRList <- c(0.4, 0.6, 0.8, 1.0)
HRList <- c(0.4, 0.6, 0.8)
colorList <- c('black', 'blue', 'darkgreen', 'red')
#labelList <- c(expression(paste(Theta, '=0.4')), expression(paste(Theta, '=0.6')), expression(paste(Theta, '=0.8')), expression(paste(Theta, '=1.0')))
labelList <- c(expression(paste(Theta, '=0.4')), expression(paste(Theta, '=0.6')), expression(paste(Theta, '=0.8')))

H12CNflux <- 0.381521
H13CNflux <- 0.403591
H12CNcont <- 0.444981
H13CNcont <- 0.442291
CO21flux <- 0.437903
CO21cont <- 0.483861
CO32flux <- 0.386186
CO32cont <- 0.442262
HCOflux <- 0.399255
HCOcont <- 0.442737

#TexGiven <- 230.0
TexGiven <- 500.0
pdf('HCN_fcov_tau.pdf', height=7)
#titleText <- expression('Solution for H'^{12}*'CN and H'^{13}*'CN flux densities')
titleText <- expression('T'[{'ex'}]*'=500 K')
contour( tauH12CN, fcov, matrix(rep(0.403591, length(fcov)*length(tauH12CN)), nrow=length(tauH12CN), ncol=length(fcov)), type='n', xlab=expression('H'^{12}*'CN optical depth'), ylab='Covering factor',  main=titleText, ylim=c(0.14, 0.7))
bestFcov <- numeric(length(HRList))
bestTau  <- numeric(length(HRList))
for( HR_index in 1:length(HRList) ){

	HR <- HRList[HR_index]
	colorName <- colorList[HR_index]
	#-------- H12CN flux
	FluxH12CN <- matrix(0, nrow=length(tauH12CN), ncol=length(fcov))
	for( tau_index in 1:length(tauH12CN) ){
		FluxH12CN[tau_index,] <- obsFluxEm(Fcov=fcov, Tau=tauH12CN[tau_index], Scont=H12CNcont, Tex=TexGiven, freq=356.734e9, r=7.4, H=HR)
	}

	#-------- H13CN flux
	FluxH13CN <- matrix(0, nrow=length(tauH13CN), ncol=length(fcov))
	for( tau_index in 1:length(tauH13CN) ){
		FluxH13CN[tau_index,] <- obsFluxEm(Fcov=fcov, Tau=tauH13CN[tau_index], Scont=H13CNcont,  Tex=TexGiven, freq=345.340e9, r=7.4, H=HR)
	}

	#-------- Draw implicit curve : Flux(fcov, tau) == observed flux
    contour( tauH12CN, fcov, FluxH12CN, add=T, levels=H12CNflux, col=colorName, labels=sprintf('%.1f', HR), lwd=2)
    contour( tauH12CN, fcov, FluxH13CN, levels=H13CNflux, add=T, col=colorName, drawlabels=F, lty=2, lwd=2)
	contour( tauH12CN, fcov, FluxH12CN, add=T, levels=c(H12CNflux-3*fluxSigmaH12CN, H12CNflux+3*fluxSigmaH12CN ), col=colorName, drawlabels=F, lwd=0.25)
	contour( tauH12CN, fcov, FluxH13CN, add=T, levels=c(H13CNflux-3*fluxSigmaH13CN, H13CNflux+3*fluxSigmaH13CN ), col=colorName, drawlabels=F, lwd=0.25)
	
	#------ most likely point
	resid <- (FluxH12CN - H12CNflux)^2 + (FluxH13CN - H13CNflux)^2
	SolutionIndex <- which(resid == min(resid), arr.ind=TRUE)
	bestTau[HR_index]  <- tauH12CN[SolutionIndex[1]]
	bestFcov[HR_index] <- fcov[SolutionIndex[2]]
	# points(bestTau[HR_index], bestFcov[HR_index], pch=20, cex=2, col=colorName)
	cat(sprintf('H/R=%.1f  Tau=%.1f  fconv=%.4f  Em=%.3e\n', HR, bestTau[HR_index], bestFcov[HR_index], FluxEm(bestFcov[HR_index], bestTau[HR_index], TexGiven, 356.734e9, 7.4, HR) ))
	
	#------ error analysis
	resid <- (FluxH12CN - H12CNflux - fluxSigmaH12CN)^2 + (FluxH13CN - H13CNflux - fluxSigmaH13CN)^2; SolutionIndex <- which(resid == min(resid), arr.ind=TRUE); cat(sprintf('H/R=%.1f  Tau=%.1f  fconv=%.4f\n', HR, tauH12CN[SolutionIndex[1]], fcov[SolutionIndex[2]] ))
	resid <- (FluxH12CN - H12CNflux + fluxSigmaH12CN)^2 + (FluxH13CN - H13CNflux - fluxSigmaH13CN)^2; SolutionIndex <- which(resid == min(resid), arr.ind=TRUE); cat(sprintf('H/R=%.1f  Tau=%.1f  fconv=%.4f\n', HR, tauH12CN[SolutionIndex[1]], fcov[SolutionIndex[2]] ))
	resid <- (FluxH12CN - H12CNflux - fluxSigmaH12CN)^2 + (FluxH13CN - H13CNflux + fluxSigmaH13CN)^2; SolutionIndex <- which(resid == min(resid), arr.ind=TRUE); cat(sprintf('H/R=%.1f  Tau=%.1f  fconv=%.4f\n', HR, tauH12CN[SolutionIndex[1]], fcov[SolutionIndex[2]] ))
	resid <- (FluxH12CN - H12CNflux + fluxSigmaH12CN)^2 + (FluxH13CN - H13CNflux + fluxSigmaH13CN)^2; SolutionIndex <- which(resid == min(resid), arr.ind=TRUE); cat(sprintf('H/R=%.1f  Tau=%.1f  fconv=%.4f\n', HR, tauH12CN[SolutionIndex[1]], fcov[SolutionIndex[2]] ))
	
	#------ Checking Flux for other lines
	cat(sprintf('CO(3-2)   flux = %.3f / %.3f\n', obsFluxEm(Fcov=bestFcov[HR_index], Tau=bestTau[HR_index], Scont=CO32cont, Tex=TexGiven, freq=345.796e9, r=7.4, H=HR), CO32flux))
	cat(sprintf('HCO+(4-3) flux = %.3f / %.3f\n', obsFluxEm(Fcov=bestFcov[HR_index], Tau=bestTau[HR_index], Scont=HCOcont, Tex=TexGiven, freq=356.734e9, r=7.4, H=HR), HCOflux))
}

for( HR_index in 1:length(HRList) ){
	points(bestTau[HR_index], bestFcov[HR_index], pch=20, cex=2, col=colorList[HR_index])
}

legend("topright", legend=labelList, col=colorList, lty=1, lwd=2)
dev.off()

tauCO32 <- bestTau* 0.1284 / 0.1508
for( HR_index in 1:length(HRList) ){
	HR <- HRList[HR_index]
	EmCO32 <- FluxEm(bestFcov[HR_index], tauCO32[HR_index],  TexGiven, 345.796e9, 7.4, HR)
	cat(sprintf('H/R=%.1f  TauCO32=%.1f Em=%.3e\n', HR, tauCO32[HR_index], EmCO32 ))
}

tauHCO <- bestTau* 0.0980 / 0.1508
for( HR_index in 1:length(HRList) ){
	HR <- HRList[HR_index]
	EmHCO <- FluxEm(bestFcov[HR_index], tauHCO[HR_index], TexGiven, 356.734e9, 7.4, HR)
	cat(sprintf('H/R=%.1f  TauHCO=%.1f Em=%.3e\n', HR, tauHCO[HR_index], EmHCO ))
}

if(0){
tauCO32 <- seq(0, 50, by=0.01)
for( HR_index in 1:length(HRList) ){
	HR <- HRList[HR_index]
	FluxCO32 <- obsFluxEm(Fcov=bestFcov[HR_index], Tau=tauCO32, Scont=CO32cont, Tex=TexGiven, freq=345.796e9, r=7.4, H=HR)
	bestTauCO32 <- tauCO32[which.min( (FluxCO32 - CO32flux)^2 )]
	cat(sprintf('H/R=%.1f  TauCO32=%.1f Em=%.3e\n', HR, bestTauCO32, FluxEm(bestFcov[HR_index], bestTauCO32, TexGiven, 345.796e9, 7.4, HR) ))
}

tauHCO <- seq(0, 50, by=0.01)
for( HR_index in 1:length(HRList) ){
	HR <- HRList[HR_index]
	FluxHCO <- obsFluxEm(Fcov=bestFcov[HR_index], Tau=tauHCO, Scont=HCOcont, Tex=TexGiven, freq=356.734e9, r=7.4, H=HR)
	bestTauHCO <- tauHCO[which.min( (FluxHCO - HCOflux)^2 )]
	cat(sprintf('H/R=%.1f  TauHCO=%.1f Em=%.3e\n', HR, bestTauHCO, FluxEm(bestFcov[HR_index], bestTauHCO, TexGiven, 356.734e9, 7.4, HR) ))
}
}