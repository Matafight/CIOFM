#Generate training set of covariates X and Z
set.seed(1)
X.tr<- matrix(rnorm(10*100),ncol = 10, nrow = 100)
Z.tr<- matrix(rnorm(15*100),ncol = 15, nrow = 100)

#Generate test set of covariates X and Z
X.te<- matrix(rnorm(10*100),ncol = 10, nrow = 100)
Z.te<- matrix(rnorm(15*100),ncol = 15, nrow = 100)
#Scale appropiately
meanX<- apply(X.tr,2,mean)
meanY<- apply(Z.tr,2,mean)
X.tr<-
Z.tr<-
X.te<-
Z.te<-
scale(X.tr, scale
scale(Z.tr, scale
scale(X.te,center
scale(Z.te,center
=
=
=
=
FALSE)
FALSE)
meanX,scale = FALSE)
meanY,scale = FALSE)
