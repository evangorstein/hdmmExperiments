#1# Loading and preparing the input data
 library(BGLR); data(mice);
 Y<-mice.pheno; X<-mice.X
 Y[is.na(Y$Biochem.Age), "Biochem.Age"] = median(Y$Biochem.Age, na.rm = TRUE) #median impute age
 y<-Y$Obesity.BMI; y<-scale(y,center=TRUE,scale=TRUE)


#2# Setting the linear predictor
 ETA<-list(  FIXED=list(~factor(GENDER)+Biochem.Age+factor(Litter),
                   data=Y,model="FIXED"),
             CAGE=list(~factor(cage),data=Y, model="BRR"),
             MRK=list(X=X,  model="BL")
       )

#3# Fitting the model
  #Move to correct folder so that files automatically created by BGLR are saved in right place
  setwd("./R/real_gwas_data/")
  set.seed(100)
  fm<-BGLR(y=y,ETA=ETA, nIter=11000, burnIn=1000)
  save(fm,file="./fm.rda")


### Box  S2 (Figure 3) #####################################################################

#1# Estimated Marker Effects & posterior SDs
  load("./fm.rda")
  bHat<- fm$ETA$MRK$b
  SD.bHat<- fm$ETA$MRK$SD.b
  plot(abs(bHat), ylab="Estimated Squared-Marker Effect",
        type="o",cex=.5,col="red",main="Marker Effects",
        xlab="Marker")
  points(abs(bHat),cex=0.5,col="blue")

#2# Predictions
  # Genomic Prediction
    gHat <- X%*%fm$ETA$MRK$b
    fixed_X = model.matrix( ~ factor(GENDER) +  Biochem.Age + factor(Litter), data=Y)[,-1]
    fixedHat <- fm$mu + fixed_X%*%fm$ETA$FIXED$b + X%*%fm$ETA$MRK$b
    yHat <- fm$yHat
    plot(fm$y~gHat,ylab="Phenotype",
         xlab="Linear predictor from genomic features", cex=0.5,
         main="Phenotype versus linear predictor from genomic features",
         xlim=range(gHat),ylim=range(fm$y));
    #abline(0, 1, col = "blue")
    plot(fm$y~fixedHat,ylab="Phenotype",
         xlab="Predicted value from all fixed effects", cex=0.5,
         main="Phenotype versus linear predictor from all fixed effects",
         xlim=range(fixedHat),ylim=range(fm$y));
    abline(0, 1, col = "blue")
    plot(fm$y~fm$yHat,ylab="Phenotype",
         xlab="Final predicted value (including random cage effects)", cex=0.5,
         main="Phenotype versus predicted value (including random cage effects)",
         xlim=range(yHat),ylim=range(fm$y));
    abline(0, 1, col = "blue")

#3# Godness of fit and related statistics
   fm$fit
   fm$varE # compare to var(y)

#4# Trace plots
  list.files()

  # Residual variance
   varE<-scan("varE.dat")
   plot(varE,type="o",col=2,cex=.5,
          ylab=expression(sigma[epsilon]^2),
          xlab="Sample",main="Residual Variance");
   abline(h=fm$varE,col=4,lwd=2);
   abline(v=fm$burnIn/fm$thin,col=4)

  # lambda (regularization parameter of the Bayesian LASSO)
   lambda<-scan("ETA_MRK_lambda.dat")
   plot(lambda,type="o",col=2,cex=.5,
          xlab="Sample",ylab=expression(lambda),
          main="Regularization parameter");
   abline(h=fm$ETA$MRK$lambda,col=4,lwd=2);
   abline(v=fm$burnIn/fm$thin,col=4)

