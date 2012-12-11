library(e1071)

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

# Select a subset with z<10 and 0.5<=E/EeV<1.5
#prSimVL <- prSim[ which(prSim$z<10. & prSim$e<1.5 & prSim$e>=0.5) , ]
#feSimVL <- feSim[ which(feSim$z<10. & feSim$e<1.5 & feSim$e>=0.5) , ]

#N_training = 100
#N_testing  = 100

#idx_pr <- sample( seq(1,length(prSimVL$e)), N_training+N_testing )
#idx_fe <- sample( seq(1,length(feSimVL$e)), N_training+N_testing )

#sampleTrain <- rbind( prSimVL[idx_pr[1:100],]   , feSimVL[idx_fe[1:100],]   )
#sampleTest  <- rbind( prSimVL[idx_pr[101:200],] , feSimVL[idx_fe[101:200],] )

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

# Straight classification without PCA
## Naive Bayes
mod.nby <- naiveBayes(comp ~ ., data=rnSimData)
prd.nby <- predict(mod.nby, rnSimData)
tbl.nby <- table(prd.nby, comp)
## SVM
#mod.svm <- svm(comp ~ ., data=rnSimData[,-8])
#prd.svm <- predict(mod.svm, rnSimData[,-8])
#tbl.svm <- table(prd.svm, comp)

# Reduction of dimensions with PCA before classifying
## PCA
pca <- princomp(rnSimData)
summary(pca) # Use the first 4 princ. componets (+94% of Var)
c1 <- t(as.matrix(pca$loadings[,1])) %*% t(rnSimData)
c2 <- t(as.matrix(pca$loadings[,2])) %*% t(rnSimData)
c3 <- t(as.matrix(pca$loadings[,3])) %*% t(rnSimData)
c4 <- t(as.matrix(pca$loadings[,4])) %*% t(rnSimData)
trnSimData <- as.data.frame( matrix( c(c1, c2, c3, c4), ncol=4, byrow=F) ) 
colnames(trnSimData) <- c('c1', 'c2', 'c3', 'c4')
## Naive Bayes
mod.pca.nby <- naiveBayes(comp ~ ., data=trnSimData)
prd.pca.nby <- predict(mod.pca.nby, trnSimData)
tbl.pca.nby <- table(prd.pca.nby, comp)
