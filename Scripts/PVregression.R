#-------- Mask by region
maskBorder <- function(x, y, startPoint, endPoint, maskSide = 'u'){
    # x          : vector of x coordinate
    # y          : vector of x coordinate
    # startPoint : 2-d vector coordinate of start point
    # endPoint   : 2-d vector coordinate of end point
    # maskSide   : u: up, d: down, l: left, r: right
    if(maskSide == 'u'){ index <- which(y > (endPoint$y - startPoint$y) / (endPoint$x - startPoint$x) * (x - startPoint$x) + startPoint$y) }
    if(maskSide == 'd'){ index <- which(y < (endPoint$y - startPoint$y) / (endPoint$x - startPoint$x) * (x - startPoint$x) + startPoint$y) }
    if(maskSide == 'r'){ index <- which(x > (endPoint$x - startPoint$x) / (endPoint$y - startPoint$y) * (y - startPoint$x) + startPoint$x) }
    if(maskSide == 'l'){ index <- which(x < (endPoint$x - startPoint$x) / (endPoint$y - startPoint$y) * (y - startPoint$x) + startPoint$x) }
    return( index )
}

#-------- projected line profile
projectLine <- function(DF, fit){
    chVeloc <- unique(DF$y)
    chNum <- length(chVeloc)
    chWidth <- mean(diff(chVeloc))
    posNum <- length(unique(DF$x))
    refVeloc <- as.vector(predict(fit, DF))[1:posNum]
    sysVeloc <- median(refVeloc)
    PVmap <- matrix(DF$z, nrow=posNum)
    cumSpec <- numeric(chNum)
    for(pos_index in 1:posNum){
        shiftCH <- round((refVeloc[pos_index] - sysVeloc) / chWidth)
        if(shiftCH >= 0){
            cumSpec[1:(chNum-shiftCH)] <- cumSpec[1:(chNum-shiftCH)] + PVmap[pos_index, (1+shiftCH):chNum]
        } else {
            cumSpec[(1-shiftCH):chNum] <- cumSpec[(1-shiftCH):chNum] + PVmap[pos_index, 1:(shiftCH+chNum)]
        }
    }
    return(data.frame(relVelocity=chVeloc-sysVeloc, flux=cumSpec))
}

#-------- Plot line with column type
plotLineColumn <- function(velocity, flux, lineColor=1, lineType=1, pointChar=20){
    lines(velocity - 0.5* median(diff(velocity)), flux, type='s', col=lineColor, lty=lineType)
    points(velocity, flux, pch=pointChar, cex=0.5, col=lineColor)
    return(0)
}

library(FITSio)
# fileName <- commandArgs(trailingOnly = T)
fileName <- 'CO32PV.fits'
setwd('.')
FITS <- readFITS(fileName)
ax1 <- axVec(1, FITS$axDat)         # Offset (arcsec)
ax2 <- axVec(2, FITS$axDat)*1.0e-3     # Velocity (m/s -> km/s)
xlab <- FITS$axDat$ctype[1]
ylab <- FITS$axDat$ctype[2]
#-------- Mask by threshold
PVimage <- FITS$imDat
threshMask <- 0.0001
pdf('CO32PV_linearRegression.pdf')
image(ax1, ax2, PVimage, xlab='Offset [arcsec]', ylab='LSR velocity [km/s]')

DF <- data.frame(
    x = rep(ax1, length(ax2)),
    y = as.vector(t(matrix(rep(ax2, length(ax1)), nrow=length(ax2)))),
    z=as.vector(PVimage))

#-------- Mask by map area
DF$z[which(DF$z < threshMask)] <- 0.0
maskPos <- c(min(DF$x), max(DF$x), 1900, 1400)
DF$z[maskBorder(DF$x, DF$y, list(x=maskPos[1], y=maskPos[3]), list(x=maskPos[2], y=maskPos[4]), 'u')] <- 0.0
lines(c(maskPos[1], maskPos[2]), c(maskPos[3],maskPos[4]))

maskPos <- c(min(DF$x), max(DF$x), 1550, 1050)
DF$z[maskBorder(DF$x, DF$y, list(x=maskPos[1], y=maskPos[3]), list(x=maskPos[2], y=maskPos[4]), 'd')] <- 0.0
lines(c(maskPos[1], maskPos[2]), c(maskPos[3],maskPos[4]))

maskPos <- c(-1.2, -1.2, 1900, 1400)
DF$z[maskBorder(DF$x, DF$y, list(x=maskPos[1], y=maskPos[3]), list(x=maskPos[2], y=maskPos[4]), 'l')] <- 0.0
lines(c(maskPos[1], maskPos[2]), c(maskPos[3],maskPos[4]))

#image(ax1, ax2, matrix(DF$z, nrow=length(ax1)), xlab='Offset [arcsec]', ylab='LSR velocity [km/s]')
fit <- lm(formula=DF$y~DF$x,  weights=DF$z)
abline(fit, lty=3, col='gray')
DF$z <- as.vector(PVimage)
DF$z[which(DF$z < -threshMask)] <- 0.0
DF$z[maskBorder(DF$x, DF$y, list(x=min(DF$x), y=1900.0), list(x=max(DF$x), y=1400.0), 'u')] <- 0.0
DF$z[maskBorder(DF$x, DF$y, list(x=min(DF$x), y=1600.0), list(x=1.0, y=1200.0), 'd')] <- 0.0
DF$z[maskBorder(DF$x, DF$y, list(x=min(-1.2), y=1900.0), list(x=max(-1.2), y=1400.0), 'l')] <- 0.0
cumSpec <- projectLine(DF, fit)
dev.off()

#-------- Plot
pdf('LineComparisonEmissoinAbsorption.pdf')
cols <- c('black', 'red', 'blue')
labels <- c('CO (J=3-2) emission', 'CO (J=2-1) absorption', 'HCN (J=4-3) absorption')
pchs <- c(20, 21, 21)
ltys <- c(1, 1, 1)
#if(0){
plot(cumSpec, pch=20, xlim=c(-420, 420), ylim=c(-0.06, 0.045), type='n', xlab='Velocity - Vsys [km/s]', ylab='Flux density [Jy]', main='Rotation-subtracted CO line profile')

GaussFit <- nls(cumSpec, formula=flux ~ a1* exp(-0.5*((relVelocity - v1)/s1)^2) + a2* exp(-0.5*((relVelocity - v2)/s2)^2) + a3* exp(-0.5*((relVelocity - v3)/s3)^2), start=list(a1=0.04, v1=-14.0, s1=59.0, a2=0.01, v2=90, s2=30, a3=0.005, v3=35, s3=10))

newDF <- data.frame(relVelocity=seq(-300, 300))
newDF$flux <- predict(GaussFit, newDF)
plotLineColumn(cumSpec$relVelocity, cumSpec$flux - predict(GaussFit, cumSpec), 'grey')

lines(newDF$relVelocity, coef(GaussFit)['a1']* exp(-0.5*((newDF$relVelocity - coef(GaussFit)['v1'])/coef(GaussFit)['s1'])^2), lty=1, col='grey')
lines(newDF$relVelocity, coef(GaussFit)['a2']* exp(-0.5*((newDF$relVelocity - coef(GaussFit)['v2'])/coef(GaussFit)['s2'])^2), lty=1, col='grey')
lines(newDF$relVelocity, coef(GaussFit)['a3']* exp(-0.5*((newDF$relVelocity - coef(GaussFit)['v3'])/coef(GaussFit)['s3'])^2), lty=1, col='grey')


abline(h=0, lty=3)
abline(v=0, lty=3)
plotLineColumn(cumSpec$relVelocity, cumSpec$flux, cols[1], ltys[1], pchs[1])

CO21 <- read.table('CO21_veloc.txt', comment.char='#')
names(CO21) <- c('veloc', 'flux')
plotLineColumn(CO21$veloc - as.numeric(coef(fit)[1]), CO21$flux - quantile(CO21$flux, 0.8), cols[2], ltys[2], pchs[2])

HCN <- read.table('HCN43_veloc.txt', comment.char='#')
names(HCN) <- c('veloc', 'flux')
plotLineColumn(HCN$veloc - as.numeric(coef(fit)[1]), HCN$flux - quantile(HCN$flux, 0.8), cols[3], ltys[3], pchs[3])

legend("bottomright", legend=labels, col=cols, pch=pchs, lty=ltys)

dev.off()