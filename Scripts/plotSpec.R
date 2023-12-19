setwd('./')
source('LineList.R')
plotLine <- function(lineSpec, lineFreq, contLevel, height=0.005, Vsys=1492.0, Vofs=0.0, Vreg=0.0, bottom=0.0){
	cvel <- 299792.458		# Velocity of light, km/s
	velPlot <- seq(-300, 300, by=100)
	# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	#cvel <- 299792.458
	

	segments( lineFreq, contLevel + height*0.2, lineFreq, contLevel + height*0.4, lwd=2 )
	segments( lineFreq, contLevel + height*0.4, lineFreq + lineSpec$offset, contLevel + height*0.7, lwd=2 )
	segments( lineFreq + lineSpec$offset, contLevel + height*0.7, lineFreq + lineSpec$offset, contLevel + height*0.9, lwd=2 )
	text( lineFreq + lineSpec$offset, contLevel + height, lineSpec$label, pos=4, offset=0.1, cex=0.8, srt=90)
	# text( freq, bottom + height, lineSpec$name, pos=1, offset=0.2, cex=1, srt=90)
	# cat(freq)
	
	#if(bottom > 0.0){
	#	freq <- (Vsys - Vofs - cvel) * lineSpec$restfreq / (Vreg + cvel)
	#	top <- bottom + height
	#	segments( min(freq), bottom, max(freq), bottom )
	#	segments( freq, rep(bottom, length(freq)), freq, rep(top, length(freq)))
	#	for( tick_index in 1:length(velPlot)){
	#		if(velPlot[tick_index] == 0){
	#			text_sd <- 'Vsys'
	#		} else {
	#			text_sd <- sprintf(' %4.0f', velPlot[tick_index])
	#		}
	#		text( freq[tick_index], top, text_sd, adj=0, cex=0.3, srt=90)
	#	}
	#	#text( mean(freq), bottom, lineSpec$name, pos=1, offset=0.2, cex=0.3)
	#}
}

plotTau <- function(DF, lineSpec, velRange, Vofs, Vreg){
	# velPlot <- seq(-100, 500, by=100)
	# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	cvel <- 299792.458
	DF$VLSR = Vofs + cvel + DF$freq/lineSpec$restfreq * (Vreg + cvel)
	plotDF <- DF[((DF$VLSR >= velRange[1]) & (DF$VLSR <= velRange[2])),]
	plotDF <- plotDF[order(plotDF$VLSR),]
	chGap <- c(0.0, diff(plotDF$VLSR))
	plot(plotDF$VLSR, plotDF$Tau, pch=20, cex=0.5, xlab='LSR Velocity [km/s]', ylab='Optical Depth', main=lineSpec$name, xlim=velRange)
	lines(plotDF$VLSR - 0.5*chGap, plotDF$Tau, type='s', lwd=0.5)
	abline(h=0, lty=2, lwd=0.5, col='grey')
	peakIndex <- which.max(plotDF$Tau)
	text_sd <- sprintf('%s ATT=%f TauMax = %f at %f km/s Flux=%f Cont=%f\n', lineSpec$name, plotDF$ATT[peakIndex], plotDF$Tau[peakIndex], plotDF$VLSR[peakIndex], plotDF$flux[peakIndex], plotDF$refSpec[peakIndex])
	cat(text_sd)
	return(plotDF)
}

freqRead <- function(prefixList){
	fileNum <- length(prefixList)
	for(file_index in 1:fileNum){
		prefix <- prefixList[file_index]
		freqFileName <- sprintf('%s.freq.txt', prefix)
		tempFreqDF  <- read.table(freqFileName);  tempFreqDF$file <- prefix
		if( file_index == 1 ){
			freqDF  <- tempFreqDF
		} else {
			freqDF <- rbind(freqDF, tempFreqDF)
		}
	}
	names(freqDF) <- c('freq', 'flux', 'file')
	return(freqDF)
}

velocRead <- function(prefixList){
	fileNum <- length(prefixList)
	spix <- numeric(fileNum)
	for(file_index in 1:fileNum){
		prefix <- prefixList[file_index]
		velocFileName <- sprintf('%s.veloc.txt', prefix)
		tempVelocDF <- read.table(velocFileName); tempVelocDF$file <- prefix
		if( file_index == 1 ){
			velocDF  <- tempVelocDF
		} else {
			velocDF <- rbind(velocDF, tempVelocDF)
		}
	}
	names(velocDF) <- c('veloc', 'flux', 'file')
	return(velocDF)
}

SpliceOverlap <- function(DF, flagDF, fileIndexA, fileIndexB, flagIndex){
	fileList <- unique(DF$file)
	chSpc <- abs(median(diff(flagDF$freq))) 	# channel spacing
	overlap_min <- min(flagDF[flagDF$fileBit == flagIndex,]$freq)
	overlap_max <- max(flagDF[flagDF$fileBit == flagIndex,]$freq)
	cat(sprintf('Overlap %.1f - %.1f GHz\n', overlap_min, overlap_max))
	fileDFa <- DF[DF$file == fileList[fileIndexA],]
	fileDFb <- DF[DF$file == fileList[fileIndexB],]
	OffsetA <- mean(fileDFa[((fileDFa$freq > overlap_min) & (fileDFa$freq < overlap_max)),]$flux)
	OffsetB <- mean(fileDFb[((fileDFb$freq > overlap_min) & (fileDFb$freq < overlap_max)),]$flux)
	DF[DF$file == fileList[fileIndexA],]$flux <- DF[DF$file == fileList[fileIndexA],]$flux - 0.5*(OffsetA - OffsetB)
	DF[DF$file == fileList[fileIndexB],]$flux <- DF[DF$file == fileList[fileIndexB],]$flux - 0.5*(OffsetB - OffsetA)
	return(DF)
}
#-------- Alignment for multiple spectral files
specAlign <- function(DF){
	fileList <- unique(DF$file)
	fileNum <- length(fileList)
	slope <- numeric(fileNum)
	contFlux <- numeric(fileNum)
	blCHList <- vector('list', fileNum)
	#-------- Check frequency overlap
	chSpc <- abs(median(diff(abs(DF$freq)))) 	# channel spacing
	freqFlag <- seq(min(DF$freq), max(DF$freq), by=chSpc)
	flagDF <- data.frame(freq=freqFlag, fileBit=rep(0, length(freqFlag)))
	for(ch_index in 1:length(DF$freq)){
		flag_index <- which( abs(flagDF$freq - DF$freq[ch_index]) < chSpc/2)
		file_index <- which(fileList == DF$file[ch_index])
		fileBit <- 2^(file_index - 1)
		flagDF$fileBit[flag_index] = flagDF$fileBit[flag_index] + fileBit
	}

	#---- check baseline (line-free) channels and get continuum spectrum
	for(file_index in 1:fileNum){
		fileDF <- DF[DF$file == fileList[file_index],]
		blFreq <- get(fileList[file_index])
		blCHList[[file_index]] <- which( ((fileDF$freq >= blFreq[1] - 0.001) & (fileDF$freq <= blFreq[2] + 0.001)) |  ((fileDF$freq >= blFreq[3] - 0.001) & (fileDF$freq <= blFreq[4] + 0.001)))
		fileDF$freq <- fileDF$freq - median(fileDF$freq)
		fit <- lm(data=fileDF[blCHList[[file_index]],], formula=flux ~ freq)
		contFlux[file_index] <- coef(fit)[1]
		slope[file_index] <- coef(fit)[2]
	}
	contFlux <- median(contFlux)
	slope <- median(slope)
	centerFreq <- mean(DF$freq)
	DF$refSpec <- contFlux + slope* (DF$freq - centerFreq)
	for(file_index in 1:fileNum){
		index <- which(DF$file == fileList[file_index])
		fileDF <- DF[DF$file == fileList[file_index],]
		scaleFact <- mean(fileDF$refSpec[blCHList[[file_index]]]) / mean(fileDF$flux[blCHList[[file_index]]])
		DF$flux[index] <- scaleFact* DF$flux[index]
	}
	#---- align overlap
	if(length(which(flagDF$fileBit == 3)) > 1){ DF <- SpliceOverlap(DF, flagDF, 1, 2, 3) }	# overlap 1 and 2 (flag=3)
	if(length(which(flagDF$fileBit == 5)) > 1){ DF <- SpliceOverlap(DF, flagDF, 1, 3, 5) }	# overlap 1 and 3 (flag=5)
	if(length(which(flagDF$fileBit == 9)) > 1){ DF <- SpliceOverlap(DF, flagDF, 1, 4, 9) }	# overlap 1 and 4 (flag=9)
	if(length(which(flagDF$fileBit == 6)) > 1){ DF <- SpliceOverlap(DF, flagDF, 2, 3, 6) }	# overlap 2 and 3 (flag=6)
	if(length(which(flagDF$fileBit == 10)) > 1){ DF <- SpliceOverlap(DF, flagDF, 2, 4, 10) }	# overlap 2 and 4 (flag=10)
	if(length(which(flagDF$fileBit == 12)) > 1){ DF <- SpliceOverlap(DF, flagDF, 3, 4, 12) }	# overlap 3 and 4 (flag=12)
	DF$ATT <- (DF$refSpec - DF$flux) / DF$refSpec
	DF$Tau <- -log(DF$flux / DF$refSpec)
	return(DF)	
}

Dopp <- function(FD, VD, fRestList){
	fileList <- unique(FD$file)
	fileNum <- length(fileList)
	D1 <- numeric(fileNum)
	D2 <- numeric(fileNum)
	for(file_index in 1:fileNum){
		index <- which(FD$file == fileList[file_index])
		fit <- lm(formula=y~x, data=data.frame(x=FD$freq[index]/fRestList[file_index], y=VD$veloc[index]))
		D1[file_index] <- coef(fit)[1]
		D2[file_index] <- coef(fit)[2]
	}
	return( matrix(c(D1, D2), ncol=2) )
}

readPlot <- function(fileList, fRestList, listofLines, plotRange, Title){
	cvel <- 299792.458		# Velocity of light, km/s
	freqDF  <- freqRead(fileList)
	velocDF <- velocRead(fileList)
	Dfacts <- Dopp(freqDF, velocDF, fRestList) - cvel
	Dofs <- median(Dfacts[,1])
	Dreg <- median(Dfacts[,2])		# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	DF <- specAlign(freqDF)
	fileNum <- length(fileList)
	#-------- Plot frame
	plot( DF$freq, DF$flux, type='n', xlim=plotRange[1:2], ylim=plotRange[3:4], xlab='Topocentric Frequency [GHz]', ylab='Flux Density [Jy]', main = Title)
	#-------- Plot spectrum
	for( file_index in 1:fileNum){
		index <- which( DF$file == fileList[file_index])
		chSep <- median(diff(DF$freq[index]))
		lines( DF$freq[index]-0.5*chSep, DF$flux[index], pch=20, col='black', type='s')	# 'steps-mid'-like lines
		points( DF$freq[index], DF$flux[index], pch=20, col='black', cex=0.2)			# data points
	}
	#-------- Plot continuum level
	lines(DF$freq, DF$refSpec, col='grey')
	
	#-------- Plot continuum level
	Vsys <- 1492.0
	for(lineID in listofLines){
		lineFreq <- (Vsys - Dofs - cvel) * lineID$restfreq / (Dreg + cvel)
		contLevel <- DF$refSpec[which.min(abs(DF$freq - lineFreq))]
		plotLine(lineID, lineFreq, contLevel, 0.005, 1492.0, Dofs, Dreg, 0.45)
		
	}
	#-------- velocity scale
	lineID <- listofLines[[1]]
	velPlot <- c(-400, -300, -200, -100, 100, 200, 300, 400)
	minFreq <- 1e9; maxFreq <- -1e9
	lineFreq <- (Vsys - Dofs - cvel) * lineID$restfreq / (Dreg + cvel)
	contLevel <- DF$refSpec[which.min(abs(DF$freq - lineFreq))]
	for( velLine in velPlot){
		lineFreq <- (Vsys + velLine - Dofs - cvel) * lineID$restfreq / (Dreg + cvel)
		minFreq <- min(lineFreq, minFreq); maxFreq <- max(lineFreq, maxFreq)
		segments( lineFreq, contLevel + 0.002, lineFreq, contLevel + 0.003, lwd=1 )
		if( abs(velLine) > 300){
			text_vel <- sprintf('%+d', velLine)
			text( lineFreq, contLevel + 0.004, text_vel, pos=4, offset=0.1, cex=0.5, srt=90)
		}
	}
	segments( minFreq, contLevel + 0.002, maxFreq, contLevel + 0.002, lwd=1 )
	
	
	for(lineID in listofLines){
		lineDF <- plotTau(DF, lineID, c(1150,2000), Dofs, Dreg)
		lineFile <- sprintf('%s.data', gsub("[ |/]", "_", lineID$name))
		write.table(lineDF[order(lineDF$VLSR),], lineFile, row.names=F, col.names=T, quote=F)
	}
	return(DF)
}

source('FileBaseline.R')
#-------- Band6 LSB
pdf('NGC1052B6LSB.pdf', width=10, height=7)
DF <- readPlot(
	# c('B6SPW0', 'B6SPW1', 'B6SPW2', 'B6SPW3'),			# fileList
	c('B6SPW2', 'B6SPW3'),
	c(217.10498, 215.220653),
	#c(230.5380, 231.90093, 217.10498, 215.220653),		# fRestList
	list(SO54, SO2_22, SO77, SO2_16, SO_34_65_54, H2S22, SiO54v0),				# listofLines
	plotRange <- c(213, 217, 0.44, 0.524),
	Title <- 'NGC 1052 Band 6 LSB'
)
dev.off()

#-------- Band6 USB
pdf('NGC1052B6USB.pdf', width=10, height=7)
DF <- readPlot(
	c('B6SPW0', 'B6SPW1'),			# fileList
	c(230.5380, 231.90093),		# fRestList
	list(CO21, H30alpha, H2O),				# listofLines
	plotRange <- c(228.4, 231.8, 0.44, 0.524),
	Title <- 'NGC 1052 Band 6 USB'
)
dev.off()

#-------- Band7 USB
pdf('NGC1052B7USB.pdf', width=10, height=7)
DF <- readPlot(
	c('B7COSPW2', 'B7COSPW3', 'B7HCNSPW0', 'B7HCNSPW1'),			# fileList
	c(356.734223, 357.921987, 354.505473, 353.622753),			# fRestList
	list(HCN43, HCN_354460, H26alpha, HCO43, HCO_356549, HCN_356256, HCO_358242, SO2_359151, SO2_357963, SO2_357926, SO2_357892, SO2_357672, SO2_357581, SO2_357388, SO2_357241, SO2_357165), # HCN_356301, HCN_356135), 				# listofLines
	plotRange <- c(351, 358, 0.38, 0.465),
	Title <- 'NGC 1052 Band 7 USB'
)
dev.off()

#-------- Band7 LSB
pdf('NGC1052B7LSB.pdf', width=10, height=7)
DF <- readPlot(
	c('B7COSPW0', 'B7COSPW1', 'B7HCNSPW2', 'B7HCNSPW3'),			# fileList
	c(345.79599, 344.916247, 342.88285, 340.71416),			# fRestList
	list(CS76, CO32, H13CN43, SO8877, H15CN43, SO_346528, SO7867, CN32a, CN32b, SO2_341276, SO2_341403, SO2_341674, SO_23_21),				# listofLines
	plotRange <- c(338, 345.2, 0.38, 0.465),
	Title <- 'NGC 1052 Band 7 LSB'
) 
dev.off()
