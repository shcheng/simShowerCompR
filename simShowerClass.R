library(e1071)
library(rpart)

# Read in data
pr.raw <- read.table('./AugerSim/simEpos_Proton_wXmax.dat')
fe.raw <- read.table('./AugerSim/simEpos_Iron_wXmax.dat')

# Keep interesting  fields
names  <- c('energy','zen','ldf.chi2','pln.chi2','rcurv','riset','s3','s4','xmax')
fields <- c(7,9,29,31,33,37,41,42,43)
pr.raw <- pr.raw[,fields]
fe.raw <- fe.raw[,fields]
colnames(pr.raw) <- names
colnames(fe.raw) <- names

# Remove NA's, 0 energy events, and 0 S4 events
idx_pr <- which( !is.na(pr.raw$x) & pr.raw$e>0 & pr.raw$ris>0 & pr.raw$rc>0)
idx_fe <- which( !is.na(fe.raw$x) & fe.raw$e>0 & fe.raw$ris>0 & fe.raw$rc>0)
pr <- pr.raw[idx_pr,]
fe <- fe.raw[idx_fe,]

# Rescale the energy to log(E/EeV) and S4 to log(S4/VEM)
## Proton
pr$energy   <- log10(pr$energy)
pr$ldf.chi2 <- log10(pr$ldf.chi2)
pr$pln.chi2 <- log10(pr$pln.chi2)
pr$rcurv    <- log10(pr$rcurv)
## Iron
fe$energy   <- log10(fe$energy)
fe$ldf.chi2 <- log10(fe$ldf.chi2)
fe$pln.chi2 <- log10(fe$pln.chi2)
fe$rcurv    <- log10(fe$rcurv)

# Add classification tag (1 for pr and 2 for fe)
prSim <- cbind( pr, as.factor(rep('pr',length(pr[,1]))) )
feSim <- cbind( fe, as.factor(rep('fe',length(fe[,1]))) )

# Properly name the data frames
colnames(prSim) <- c(names, 'comp')
colnames(feSim) <- c(names, 'comp')

simData <- rbind(prSim, feSim)


ztab <- matrix(c(22, 29, 43, 47, 55, 58.75), ncol=2, byrow=T) 
etab <- matrix(c(18.0, 18.66, 18.66, 19.15, 19.15, 19.65), ncol=2, byrow=T)

pdf("xmax_s4.pdf")
par(mfrow=c(3,3))
for( i in 1:(dim(ztab)[1]) ) { 
  for( j in 1:(dim(etab)[1]) ) {
    db <- simData[which(simData$z>=ztab[i,1] & simData$z<ztab[i,2] & 
                        simData$e>=etab[j,1] & simData$e<etab[j,2]),]
    zmin <- as.character(ztab[i,1])
    zmax <- as.character(ztab[i,2])
    emin <- as.character(etab[j,1])
    emax <- as.character(etab[j,2])
    title <- paste('Z [', zmin, ',', zmax, ') logE [', emin, ',', emax, ')', sep="")
    xlabel <- 'Xmax [g/cm^2]'
    ylabel <- 'S4 [VEM]'
    plot(db$xmax, db$s4, pch=21, ,cex=0.5, col=db$comp, main=title, xlab=xlabel, ylab=ylabel)
  }
}
dev.off()

pdf("riset_s4.pdf")
par(mfrow=c(3,3))
for( i in 1:(dim(ztab)[1]) ) { 
  for( j in 1:(dim(etab)[1]) ) {
    db <- simData[which(simData$z>=ztab[i,1] & simData$z<ztab[i,2] & 
                        simData$e>=etab[j,1] & simData$e<etab[j,2]),]
    zmin <- as.character(ztab[i,1])
    zmax <- as.character(ztab[i,2])
    emin <- as.character(etab[j,1])
    emax <- as.character(etab[j,2])
    title <- paste('Z [', zmin, ',', zmax, ') logE [', emin, ',', emax, ')', sep="")
    xlabel <- 'risetime [nsec]'
    ylabel <- 'S4 [VEM]'
    plot(db$riset, db$s4, pch=21, ,cex=0.5, col=db$comp, main=title, xlab=xlabel, ylab=ylabel)
  }
}
dev.off()

pdf("LDFchi2_s4.pdf")
par(mfrow=c(3,3))
for( i in 1:(dim(ztab)[1]) ) { 
  for( j in 1:(dim(etab)[1]) ) {
    db <- simData[which(simData$z>=ztab[i,1] & simData$z<ztab[i,2] & 
                        simData$e>=etab[j,1] & simData$e<etab[j,2]),]
    zmin <- as.character(ztab[i,1])
    zmax <- as.character(ztab[i,2])
    emin <- as.character(etab[j,1])
    emax <- as.character(etab[j,2])
    title <- paste('Z [', zmin, ',', zmax, ') logE [', emin, ',', emax, ')', sep="")
    xlabel <- 'log( LDF chi^2 ) [a.u.]'
    ylabel <- 'S4 [VEM]'
    plot(db$ldf.chi2, db$s4, pch=21, ,cex=0.5, col=db$comp, main=title, xlab=xlabel, ylab=ylabel)
  }
}
dev.off()

pdf("PLANEchi2_s4.pdf")
par(mfrow=c(3,3))
for( i in 1:(dim(ztab)[1]) ) { 
  for( j in 1:(dim(etab)[1]) ) {
    db <- simData[which(simData$z>=ztab[i,1] & simData$z<ztab[i,2] & 
                        simData$e>=etab[j,1] & simData$e<etab[j,2]),]
    zmin <- as.character(ztab[i,1])
    zmax <- as.character(ztab[i,2])
    emin <- as.character(etab[j,1])
    emax <- as.character(etab[j,2])
    title <- paste('Z [', zmin, ',', zmax, ') logE [', emin, ',', emax, ')', sep="")
    xlabel <- 'log( Plane chi^2 ) [a.u.]'
    ylabel <- 'S4 [VEM]'
    plot(db$pln.chi2, db$s4, pch=21, ,cex=0.5, col=db$comp, main=title, xlab=xlabel, ylab=ylabel)
  }
}
dev.off()

pdf("rcurv_s4.pdf")
par(mfrow=c(3,3))
for( i in 1:(dim(ztab)[1]) ) { 
  for( j in 1:(dim(etab)[1]) ) {
    db <- simData[which(simData$z>=ztab[i,1] & simData$z<ztab[i,2] & 
                        simData$e>=etab[j,1] & simData$e<etab[j,2]),]
    zmin <- as.character(ztab[i,1])
    zmax <- as.character(ztab[i,2])
    emin <- as.character(etab[j,1])
    emax <- as.character(etab[j,2])
    title <- paste('Z [', zmin, ',', zmax, ') logE [', emin, ',', emax, ')', sep="")
    xlabel <- 'log( Rcurv ) [a.u.]'
    ylabel <- 'S4 [VEM]'
    plot(db$rcurv, db$s4, pch=21, ,cex=0.5, col=db$comp, main=title, xlab=xlabel, ylab=ylabel)
  }
}
dev.off()

# Normalize simData
colNorm <- function( c ) {
    return ((c - mean(c))/sd(c))
}

rnSimData <- as.data.frame( apply(simData[,c(1:6,8)], 2, colNorm) )
comp <- simData[,10]

etab <- c(17.4, 18.10, 18.65, 19.15, 19.65, 20.10, 20.65, 21.20)
emid <- NULL
for(i in 1:(length(etab)-1)) {
  emid <- c(emid, (etab[i]+etab[i+1])/2)
}
classErr.nby <- NULL
prClassErr.nby <- NULL
feClassErr.nby <- NULL
classErr.svm <- NULL
prClassErr.svm <- NULL
feClassErr.svm <- NULL
classErr.crt <- NULL
for(i in 1:(length(etab)-1)) {
  elow_char <- as.character(etab[i])
  ehig_char <- as.character(etab[i+1])
  print(">>>")
  print(paste(elow_char, ehig_char, sep=" "))
  iCov <- which(simData$energy>=etab[i] & simData$energy<etab[i+1])
  icomp <- comp[iCov] 
  irnData <- rnSimData[iCov,]
  pca <- princomp(irnData)
  c1 <- t( t(as.matrix(pca$loadings[,1])) %*% t(irnData) )
  c2 <- t( t(as.matrix(pca$loadings[,2])) %*% t(irnData) )
  c3 <- t( t(as.matrix(pca$loadings[,3])) %*% t(irnData) )
  c4 <- t( t(as.matrix(pca$loadings[,4])) %*% t(irnData) )
  pcData <- as.data.frame( matrix(c(c1, c2, c3, c4), ncol=4, byrow=F) )
  colnames(pcData) <- c('comp1', 'comp2', 'comp3', 'comp4')
  ## Create training set and testing set
  rndIdx <- sample.int(length(icomp))
  shuffledData   <- pcData[rndIdx,]
  shuffledTarget <- icomp[rndIdx]
  critIdx <- floor(length(shuffledData[,1])*0.5)
  trainData <- shuffledData[1:critIdx,]
  trainTarget <- shuffledTarget[1:critIdx]
  testData <- shuffledData[(critIdx+1):length(shuffledData[,1]),]
  testTarget <- shuffledTarget[(critIdx+1):length(shuffledData[,1])]
  ## Naive Bayes Classifier
  nby.mod <- naiveBayes(trainTarget ~ ., data=trainData)
  nby.prd <- predict(nby.mod, testData)
  nby.tab <- table(testTarget, nby.prd)
  print(nby.tab)
  classErr.nby <- c(classErr.nby, 1-sum(nby.prd==testTarget)/length(nby.prd))
  prClassErr.nby <- c(prClassErr.nby, nby.tab[1,2]/sum(nby.tab[1,]))
  feClassErr.nby <- c(feClassErr.nby, nby.tab[2,1]/sum(nby.tab[2,]))
  ## SVM
  svm.mod <- svm(trainTarget ~ ., data=trainData)
  svm.prd <- predict(svm.mod, testData)
  svm.tab <- table(testTarget, svm.prd)
  print(svm.tab)
  classErr.svm <- c(classErr.svm, 1-sum(svm.prd==testTarget)/length(svm.prd))
  prClassErr.svm <- c(prClassErr.svm, svm.tab[1,2]/sum(svm.tab[1,]))
  feClassErr.svm <- c(feClassErr.svm, svm.tab[2,1]/sum(svm.tab[2,]))
  ## rpart
  crt.mod <- rpart(trainTarget ~ ., data=trainData)
  crt.reg <- predict(crt.mod, testData)
  crt.prd <- NULL
  crt.prd[(crt.reg[,1]>=crt.reg[,2])] = 'pr'
  crt.prd[(crt.reg[,1]<crt.reg[,2])]  = 'fe'
  crt.prd <- as.factor(crt.prd)
  #crt.tab <- table(testTarget, crt.prd)
  classErr.crt <- c(classErr.crt, 1-sum(crt.prd==testTarget)/length(crt.prd))
}

# Plot the classification errors for each algorithm
## plot a empty frame
plot(emid, classErr.nby, xlim=c(floor(etab[1]), 21.5), 
     ylim=c(0.1, 0.4), xlab='log(E)', ylab='Classification Error',
     main="Error in Classification")
## Add line markers for the energy ranges
for( i in 1:length(etab) ) {
  abline(v=etab[i], col='green', lty=2)
  energy.char <- as.character(etab[i])
  text(etab[i], 0.1, energy.char)
}
## Add horizontal eye guides 
abline(h=seq(0.1,0.4,0.05), col='gray', lty=1)
## Add the actual plots
points(emid, classErr.nby, ty='b', pch=1, col='black')
points(emid, classErr.crt, ty='b', pch=2, col='red')
points(emid, classErr.svm, ty='b', pch=5, col='blue')
legend('topright', c("Naive Bayes", "CART", "SVM"), pch=c(1,2,5), 
       col=c('black','red','blue'), lty=rep(1,3), bg='white')
