setwd('./')
source('LineList.R')
plotLine <- function(lineSpec, velPlot, Vofs, Vreg, bottom, height){
	# velPlot <- seq(-100, 500, by=100)
	# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	cvel <- 299792.458
	top <- bottom + height
	freq <- (velPlot - Vofs - cvel) * lineSpec$restfreq / (Vreg + cvel)
	cat(freq)
	segments( min(freq), bottom, max(freq), bottom )
	segments( freq, rep(bottom, length(freq)), freq, rep(top, length(freq)))
	for( tick_index in 1:length(velPlot)){
		if(velPlot[tick_index] == 0){
			text_sd <- 'Vsys'
		} else {
			text_sd <- sprintf(' %4.0f', velPlot[tick_index])
		}
		text( freq[tick_index], top, text_sd, adj=0, cex=0.3, srt=90)
	}
	text( mean(freq), bottom, lineSpec$name, pos=1, offset=0.2, cex=0.3)
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

specAlign <- function(DF){
	fileList <- unique(DF$file)
	fileNum <- length(fileList)
	slope <- numeric(fileNum)
	contFlux <- numeric(fileNum)
	blCHList <- vector('list', fileNum)
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

readPlot <- function(fileList, fRestList, listofLines, labelY, plotRange, Title){
	cvel <- 299792.458		# Velocity of light, km/s
	freqDF  <- freqRead(fileList)
	velocDF <- velocRead(fileList)
	Dfacts <- Dopp(freqDF, velocDF, fRestList) - cvel
	Dofs <- median(Dfacts[,1])
	Dreg <- median(Dfacts[,2])		# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	DF <- specAlign(freqDF)
	fileNum <- length(fileList)
	plot( DF$freq, DF$flux, type='n', xlim=plotRange[1:2], ylim=plotRange[3:4], xlab='Frequency [GHz]', ylab='Flux Density [Jy]', main = Title)
	for( file_index in 1:fileNum){
		index <- which( DF$file == fileList[file_index])
		chSep <- median(diff(DF$freq[index]))
		lines( DF$freq[index]-0.5*chSep, DF$flux[index], pch=20, col=file_index, type='s')
		points( DF$freq[index], DF$flux[index], pch=20, col='black', cex=0.1)
	}
	lines(DF$freq, DF$refSpec, col='grey')
	velPlot <- c(1000, 1500, 2000)
	index <- 1
	for(lineID in listofLines){
		plotLine(lineID, velPlot, Dofs, Dreg, labelY[index], 0.001)
		index <- index + 1
	}
	for(lineID in listofLines){
		lineDF <- plotTau(DF, lineID, c(1000,2000), Dofs, Dreg)
		lineFile <- sprintf('%s.data', gsub("[ |/]", "_", lineID$name))
		write.table(lineDF[order(lineDF$VLSR),], lineFile, row.names=F, col.names=T, quote=F)
	}
	return(DF)
}

source('FileBaseline.R')

pdf('NGC1052B6.pdf')
DF <- readPlot(
	c('B6SPW0', 'B6SPW1', 'B6SPW2', 'B6SPW3'),			# fileList
	c(230.5380, 231.90093, 217.10498, 215.220653),			# fRestList
	list(CO21, H2O, SO77, SO2_16, SO54, SO_34_65_54, SO2_22),				# listofLines
	labelY <- c(0.49, 0.488, 0.505, 0.517, 0.512, 0.507, 0.502),
	plotRange <- c(213, 231.8, 0.44, 0.52),
	Title <- 'NGC 1052 Band 6'
)
dev.off()

pdf('NGC1052B7.pdf')
#-------- Band7 USB
DF <- readPlot(
	c('B7COSPW2', 'B7COSPW3', 'B7HCNSPW0', 'B7HCNSPW1'),			# fileList
	c(356.734223, 357.921987, 354.505473, 353.622753),			# fRestList
	list(HCO43, HCN43,HCO_356549, HCN_356256),				# listofLines
	labelY <- c(0.445, 0.445, 0.45, 0.455),
	plotRange <- c(351, 358, 0.38, 0.46),
	Title <- 'NGC 1052 Band 7 USB'
)

#-------- Band7 LSB
DF <- readPlot(
	c('B7COSPW0', 'B7COSPW1', 'B7HCNSPW2', 'B7HCNSPW3'),			# fileList
	c(345.79599, 344.916247, 342.88285, 340.71416),			# fRestList
	list(CO32, H13CN43, SO8877, SO_346528, CS76, SO78676, CN32a, CN32b, HCS87, NaCN),				# listofLines
	labelY <- c(0.445, 0.45, 0.445, 0.45, 0.447, 0.447, 0.457, 0.452, 0.452, 0.447),
	plotRange <- c(338, 345.2, 0.38, 0.46),
	Title <- 'NGC 1052 Band 7 LSB'
) 
dev.off()
