#-------- Plot line with column type
plotLineColumn <- function(velocity, flux, lineColor=1, lineType=1, pointChar=20){
    lines(velocity - 0.5* median(diff(velocity)), flux, type='s', col=lineColor, lty=lineType)
    points(velocity, flux, pch=pointChar, cex=0.5, col=lineColor)
    return(0)
}

# fileName <- commandArgs(trailingOnly = T)
prefixList <- c('-0.44_-2.44', '-0.72_-0.12', '-1.00_-1.00', '-1.52_-4.56', '-2.76_-0.76', '-2.76_0.12', '-3.72_-0.64', '-3.84_-1.16', '-4.16_-3.60', '-5.08_-4.48', '-5.68_-3.68', '1.36_1.68', '1.52_-0.32', '2.56_1.16', '2.60_1.48', '2.80_0.52', '2.80_1.80', '3.64_1.52', '4.16_-0.40', '6.64_-0.68', '7.00_0.24', '7.88_1.76', '9.20_3.72', 'CND_integFlux')
lineRangeList <- list(
c(1420, 1540),  # -0.44_-2.44
c(1890, 2000),  # -0.72_-0.12
c(1470, 1760),  # -1.00_-1.00
c(1650, 1760),  # -1.52_-4.56
c(1500, 1660),  # -2.76_-0.76
c(1550, 1620),  # -2.76_0.12
c(1980, 2030),  # -3.72_-0.64
c(2080, 2160),  # -3.84_-1.16
c(1630, 1680),  # -4.16_-3.60
c(2070, 2130),  # -5.08_-4.48
c(2020, 2090),  # -5.68_-3.68
c(1160, 1240),  # 1.36_1.68
c(1570, 1620),  # 1.52_-0.32
c(1250, 1360),  # 2.56_1.16
c(1280, 1350),  # 2.60_1.48
c(1250, 1340),  # 2.80_0.52
c(1220, 1450),  # 2.80_1.80
c(1210, 1460),  # 3.64_1.52
c(1240, 1400),  # 4.16_-0.40
c(1930, 1970),  # 6.64_-0.68
c(1400, 1500),  # 7.00_0.24
c(2030, 2080),  # 7.88_1.76
c(2010, 2180),  # 9.20_3.72
c(1190, 1660))  # CND_integFlux
setwd('.')
fileNum <- length(prefixList)

setLineDF <- function(fileName, lineRange){
    DF <- read.table(fileName, comment.char='#')
    names(DF) <- c('veloc', 'flux')
    lineIndex <- which((DF$veloc > lineRange[1]) & (DF$veloc < lineRange[2]))
    freeIndex <- which((DF$veloc < lineRange[1]) | (DF$veloc > lineRange[2]))
    specOffset <- median(DF$flux[freeIndex])
    DF$flux <- DF$flux - specOffset
    return(list(DF, lineIndex, freeIndex))
}

plotLineDF <- function(DF, fileName, lineIndex, freeIndex, colorLine){
    plotLineColumn(DF$veloc, DF$flux - median(DF$flux), colorLine)
    chWidth <- median(abs(diff(DF$veloc)))
    sumFlux <- sum(DF$flux[lineIndex])
    errFlux <- sqrt(length(lineIndex))* sd(DF$flux[freeIndex])
    peakVel <- DF$veloc[lineIndex[which.max(DF$flux[lineIndex])]]
    SNR <- sumFlux/errFlux
    text_sd <- sprintf('%s : Peak = %f km/s, Integ. Flux = %f (%f) Jy km/s SNR=%f\n', fileName, peakVel, chWidth* sumFlux, chWidth*errFlux, SNR)
    cat(text_sd)
    return(0)
}

cols <- c('red', 'blue')

for(file_index in 1:fileNum){
    file21 <- sprintf('CO21_%s.txt', prefixList[file_index]) 
    file32 <- sprintf('CO32_%s.txt', prefixList[file_index]) 
    pdf(sprintf('CO%s.pdf', prefixList[file_index]))
    lineRange <- lineRangeList[[file_index]]
    setList <- setLineDF(file21, lineRange); DF21 <- setList[[1]]; lineIndex21 <- setList[[2]]; freeIndex21 <- setList[[3]]
    setList <- setLineDF(file32, lineRange); DF32 <- setList[[1]]; lineIndex32 <- setList[[2]]; freeIndex32 <- setList[[3]]
    plot(DF21, pch=20, xlim=lineRange+c(-250,250), ylim=c(-0.1*max(DF32$flux), 1.2*max(DF32$flux)), type='n', xlab='LSR Velocity [km/s]', ylab='Flux density [Jy]', main='CO line profiles')
    abline(h=0, col='grey')
    abline(v=lineRange[1], col='grey')
    abline(v=lineRange[2], col='grey')
    plotLineDF(DF21, file21, lineIndex21, freeIndex21, cols[1])
    plotLineDF(DF32, file32, lineIndex32, freeIndex32, cols[2])
    dev.off()
}
