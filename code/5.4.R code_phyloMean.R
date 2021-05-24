gls.Ymean<-function(Y,Sigma){
      n<-length(Y)
#Standardize sigma to reduce rounding errors
#      tr<-sum(diag(Sigma))
#      Sigma<-n*Sigma/tr
#Input
      invSigma<-solve(Sigma)
#pgls mean and variance
      X1<-rep(1,n)
      q<-2          # correct if multivariate!!
      C1<-solve(t(X1)%*%invSigma%*%X1)
      Y_PGLSmean<-C1%*%t(X1)%*%invSigma%*%Y
      Y_PGLSdeviations = Y - c(Y_PGLSmean)
      Y_PGLSvariance = (t(Y_PGLSdeviations)%*%invSigma%*%Y_PGLSdeviations)/(n-1)
      SE_Y_mean = sqrt(Y_PGLSvariance/n)
#Save model
      results<-cbind(Y_PGLSmean,SE_Y_mean,Y_PGLSvariance)
      colnames(results)<-c("Ymean","YSE","Y_PGLSvariance")
return(results)
}
