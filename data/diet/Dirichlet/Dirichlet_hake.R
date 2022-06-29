library(VGAM)   #fits Dirichlet Dist.
library(crayon) #color output to console

path <- "data/diet/Dirichlet/plots/"  # path for all plots

x <- read.csv("data/diet/Dirichlet/hake_for_dirichlet.csv")

#read data   
npredspp  <- nlevels(as.factor(x$Predator))
nprey     <- ncol(x) - 2
predatornames <- unique(x$Predator)

#--------------CHECK ME BEFORE RUNNING-----------------
# What graphs do you want to produce? Chose whatplot 1-4
whatplot <- 1  # (fitted marginal betas for all predators)
# whatplot <- 2  # (bootstrapped data histogram)
# whatplot <- 3  # (A larger plot of marginal beta dist for a particular predator) 
# whatplot <- 4  # (shows marginal beta, bootstrap histogram and simple average for a particular predator-prey interaction)

FirstPredToDo <- 1
LastPredToDo <- as.numeric(npredspp)   
whatplot3pred <- as.numeric(npredspp) 
whatplot4pred <-  as.numeric(npredspp)  # currently set for 15th predator in file
whatplot4prey <- 1  # set for 1st predator in file

# Plot options
options(digits = 20)  
nfishavg <- 10  # number of stomachs to average
nbootstraps <- 1000

if (whatplot == 3) {
  FirstPredToDo = whatplot3pred
  LastPredToDo = whatplot3pred
}

if (whatplot == 4) {
  FirstPredToDo = whatplot4pred
  LastPredToDo = whatplot4pred
}


#-------------------------------INITIALIZATION-----------------------------------
#dirichlet fit will fail if you feed it 0 or 1, so squish the values a bit on both sides
dat2 <- cbind(x[,1:2],(0.999*(0.00001 + x[,3:(nprey+2)]/rowSums(x[,3:(nprey+2)]))))

#set up matrices to hold results
mat_names   <- list(c(colnames(dat2[,3:(nprey+2)])),c(levels(dat2$Predator)))
raw_wt_avg  <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
mle         <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
betavar     <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
lower95     <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
upper95     <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
mode_       <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
beta_a      <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
beta_b      <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
simpleaverage = matrix(nrow=npredspp,ncol=nprey)
the_mode = matrix(nrow=nprey,ncol=npredspp)
bootaverage = matrix(nrow=nprey,ncol=npredspp)
bootmode = matrix(nrow=nprey,ncol=npredspp)


#-------------------------------BOOTSTRAPPING-----------------------------------
for(i in FirstPredToDo:LastPredToDo){         
  
  #dat3 <- subset(dat2[dat2$Predator==levels(dat2$Predator)[i],]) #dat3 contains data for only pred i
  
  #thispred = levels(as.factor(dat2$Predator))[i]
  thispred = predatornames[i]
  dat3 = subset(dat2[dat2$Predator==thispred,])
  
  simpleaverage[i,]=colMeans(dat3[3:(nprey+2)])                  #verified correct
  
  nstomachs   <- nrow(dat3)
  stor        <- matrix(nrow=nbootstraps,ncol=nprey)
  tempstor <- matrix(nrow=nbootstraps,ncol=nprey)
  if (nstomachs<nfishavg) {print(paste("Warning: there are only ", nstomachs, " stomachs for ",predatornames[i]))}
  
  for(j in 1:nbootstraps){
    draw      <- dat3[sample(nstomachs,nfishavg,replace=T),]
    
    #-----
    #tempstor[j,] <- apply(draw[,3:(nprey+2)],2,weighted.mean,w=draw$wgt) 
    drawsum =  colSums(draw[,3:(nprey+2)]) / length(draw[,1])
    tempstor[j,]=drawsum
    
    stor[j,] = tempstor[j,] / sum(tempstor[j,])
    
  }
  
  for (k in 1:nprey){
    bootaverage[k,i]=mean(stor[,k])
    bootmode[k,i]=as.numeric(names(sort(-table(stor[,k])))[1])
  }
  #---------------------------FIT TO DIRICHLET DIST-------------------------------
  
  cat(paste(bgGreen$white$bold(as.character(predatornames[i])),"\n"))
  print("trying to fit")
  fit <- vglm(stor~1,dirichlet)
  print("fit worked")
  k   <- Coef(fit)
  
  #obtain marginal beta distribution parameters
  betaEst <- matrix(nrow=nprey,ncol=2)
  for(m in 1:nprey){
    betaEst[m,1] <- k[m]
    betaEst[m,2] <- sum(k)-k[m]
  }
  
  #find pdf of beta for each prey item looping over a range of diet proportions
  LogDietCategories = 10000
  prop = exp(seq(log(0.0000000001), log(0.999), length.out = LogDietCategories)) 
  
  bdist   <- data.frame(matrix(nrow=length(prop),ncol=nprey))
  
  print("trying to populate bdist")
  for(n in 1:length(prop)){      
    bdist[n,] <- dbeta(prop[n],betaEst[,1],betaEst[,2])
  }
  print("population successful")
  
  #-----------------------COMPUTE AND STORE THE STATISTICS----------------------------------
  #mean of raw data
  #raw_wt_avg[,i] <- apply(dat3[,3:(nprey+2)],2,weighted.mean,w=dat3$wgt)
  drawsum2 = dat3[,3:(nprey+2)]
  drawsum2 =  colSums(dat3[,3:(nprey+2)]) / length(dat3[,1])
  raw_wt_avg[,i]=drawsum2
  
  #maximum likelihood estimate
  Y <- apply(bdist,2,max)
  I <- apply(bdist,2,which.max)
  maxes <- prop[c(I)]
  mle[,i] <- maxes   
  
  #variance
  betavar[,i] <- betaEst[,1]*betaEst[,2]/((betaEst[,1]+betaEst[,2])^2 * (betaEst[,1]+betaEst[,2]+1))
  
  #95% confidence intervals
  for(p in 1:nprey){
    lower95[p,i] <- qbeta(.025,betaEst[p,1],betaEst[p,2])
    upper95[p,i] <- qbeta(.975,betaEst[p,1],betaEst[p,2])
  }
  
  #Beta Parameters
  beta_a[,i] <- t(betaEst)[1,]
  beta_b[,i] <- t(betaEst)[2,] 
  
  #--------------------------------PLOTS------------------------------------------
  #identify which prey were actually eaten
  preyeaten=which(colSums(tempstor)>0)
  
  #generate frequency data sets
  bins <- seq(0,1.1,.01)  
  binmids <- seq(.005,1.1,.01)   
  N <- data.frame(matrix(nrow=length(bins)-1,ncol=nprey))
  N2 <- data.frame(matrix(nrow=length(bins)-1,ncol=nprey))
  for(r in 1:nprey){
    N[,r]   <- (hist(stor[,r],breaks=bins,plot=FALSE)$density)/100
    N2[,r]  <- (hist(dat3[,r+2],breaks=bins,plot=FALSE)$density)/100
  }
  
  #................................whatplot1..........................................
  if(whatplot==1){
    
    if(i%%10==1|i==FirstPredToDo){
      
      FigName=paste0(path, "MargBeta1.tif")
      if(i>10){FigName=paste0(path, "MargBeta2.tif")}
      if(i>20){FigName=paste0(path, "MargBeta3.tif")}
      if(i>30){FigName=paste0(path, "MargBeta4.tif")}
      if(i>40){FigName=paste0(path, "MargBeta5.tif")}
      if(i>50){FigName=paste0(path, "MargBeta6.tif")}
      if(i>60){FigName=paste0(path, "MargBeta7.tif")}
      
      tiff(FigName,width=7, height=7, units="in", res=600, compression = "lzw")
      
      par(mfcol=c(5,2))
      par(mar=c(2,2,2,2))
      par(omi=c(.5,.5,0,.5))
    }
    
    #scale the y axis according to the maximum beta that ISNT in the leftmost category
    bigpreyeaten=preyeaten[I[preyeaten]!=1]
    
    #only show the prey whose mode is greater than the smallest category
    plot(prop,bdist[,1],type='n',col="Black",ylim=c(0,max(bdist[,bigpreyeaten])),xlim=c(0,1),main=paste(dat3$Predator[1]," (prey items=",length(bigpreyeaten),")",sep=""))
    for(y in 1:nprey){
      
      printprey=FALSE
      if (sum(tempstor[,y])>0){printprey=TRUE}
      
      if(printprey==TRUE){
        lines(prop,bdist[,y],col="Black")
      }
    }
    #legend(x=.75,y=(max(bdist)*.8),c(colnames(dat3[,3:(nprey+2)])),col=c(1,2,3),pch=16)
    
    if(i%%10==0|i==LastPredToDo){   
      mtext("Contribution to diet",1,outer=T,adj=.5,1.5)   
      mtext("Likelihood",2,outer=T,adj=.5,1.5)
      dev.off()
    }
  }
  
  #................................whatplot2..........................................
  if(whatplot==2){
    
    if(i%%10==1|i==FirstPredToDo){
      FigName=paste0(path, "Bootstrap1.tif")
      if(i>10){FigName=paste0(path, "Bootstrap2.tif")}
      if(i>20){FigName=paste0(path, "Bootstrap3.tif")}
      if(i>30){FigName=paste0(path, "Bootstrap4.tif")}
      if(i>40){FigName=paste0(path, "Bootstrap5.tif")}
      if(i>50){FigName=paste0(path, "Bootstrap6.tif")}
      if(i>60){FigName=paste0(path, "Bootstrap7.tif")}
      
      tiff(FigName,width=7, height=7, units="in", res=600, compression = "lzw")
      
      par(mfcol=c(5,2))
      par(mar=c(2,2,2,2))
      par(omi=c(.5,.5,0,.5))
    }
    
    plot(binmids[N[,1]>0],N[N[,1]>0,1],ylim=c(0,.4),xlim=c(0,1),type='n',main=dat3$Predator[1])
    for(xx in 1:nprey){
      lines(binmids[N[,xx]>0],N[N[,xx]>0,xx],col=xx)
      #points(binmids[N2[,xx]>0],N2[N2[,xx]>0,xx],col=xx,pch=16)
      lines(c(raw_wt_avg[xx,i],raw_wt_avg[xx,i]),c(0,.4),col=xx,lty=2)  
    }
    
    if(i%%10==0|i==LastPredToDo){    
      mtext("Contribution to diet",1,outer=T,adj=.5,1.5)   
      mtext("Frequency",2,outer=T,adj=.5,1.5)
      dev.off()
    }
  }
  
  #................................whatplot4..........................................
  if(whatplot==4){
    
    if(i==whatplot4pred){
      preyname = gsub("_","",names(dat3)[2+whatplot4prey])
      FigName = paste(path, gsub(" ","",predatornames[i]),"_",preyname,".tif",sep="")
      tiff(FigName,width=7, height=7, units="in", res=600, compression = "lzw")
      
      par(mar=c(2,2,2,2))
      par(omi=c(.5,.5,0,.5))
      
      bigpreyeaten=preyeaten[I[preyeaten]!=1]
      
      #firsttime = TRUE
      
      y=whatplot4prey
      lowdone=0
      highdone=0
      BdistForPlotting=bdist[,y]/sum(bdist[,y])
      bdistsum = 0
      lastBdist = 0
      ii=0
      #determine 99 percentile for x limit of plot
      #could also modify this to plot 95% CI if you were so inclined
      while (ii < LogDietCategories) {
        ii = ii + 1
        bdistsum = bdistsum + BdistForPlotting[ii]
        
        if (lowdone == FALSE) {
          if (bdistsum > (0.001 * sum(BdistForPlotting))){
            X01 = prop[ii]
            print(ii)
            lowdone=TRUE
          }
        }
        
        if (BdistForPlotting[ii]>lastBdist){
          lastBdist = BdistForPlotting[ii]
          X50 = prop[ii]
        }
        
        if (highdone == FALSE){
          if (bdistsum > (0.999 * sum(BdistForPlotting))){
            X99 = prop[ii]
            print(ii)
            highdone=TRUE
          }
        }
      }     
      
      X01=X01*.9
      X99=X99+(0.1*X01)
      
      plot(prop,BdistForPlotting,col="Black",type="l",ylim=c(0,max(BdistForPlotting)),xlab="",ylab="",xlim=c(X01,X99),main=gsub(" ","",paste(dat3$Predator[1],preyname,sep="_")))
      mtext("Probability",side=2,outer=FALSE,line=2.5)
      mtext("Contribution to diet",side=1,outer=FALSE,line=2.5)
      
      #histogram
      a=hist(stor[,y],breaks=50,plot=FALSE)     
      par(new=TRUE)
      plot(a,xlim=c(X01,X99),axes=FALSE,col="light gray",border="white",xlab="",ylab="",main="")
      axis(side=4)
      mtext("Count",side=4,outer=FALSE,line=2.1)
      
      par(new=TRUE)
      plot(prop,BdistForPlotting,col="Black",type="l",ylim=c(0,max(BdistForPlotting)),axes=FALSE,xlim=c(X01,X99))
      
      firsttime = FALSE
      
      #vertical lines for simple average
      lines(c(bootaverage[y,i],bootaverage[y,i]),c(0,2000),col="Black",lty=2)
      
      #vertical line for mode
      lines(c(X50,X50),c(0,2000),col="Black",lty=3)
      
      dev.off() 
    }
  }
  
  #................................whatplot3.......................................... 
  if(whatplot==3){ 
    if(i == whatplot3pred){
      FigName=paste(path, dat3$Predator[1], "Beta.tif")
      FigName=gsub(" ","",FigName)
      tiff(FigName,width=7, height=7, units="in", res=600, compression = "lzw")
      
      par(mar=c(2,2,2,2))
      par(omi=c(.5,.5,0,.5))
      
      #scale the y axis according to the maximum beta that ISNT in the leftmost category
      bigpreyeaten=preyeaten[I[preyeaten]!=1] 
      
      firsttime = TRUE
      
      for(y in 1:nprey){
        
        printprey=FALSE
        if (mle[y,i]>0.00000001){printprey=TRUE}
        #BdistForPlotting=bdist[,y]/sum(bdist[,y])
        BdistForPlotting=bdist[,y]
        if(printprey==TRUE){
          
          if(firsttime==TRUE){
            plot(prop,BdistForPlotting,type="n",col="Black",ylim=c(0,max(bdist[bigpreyeaten])),xlim=c(0,0.5),main=paste(dat3$Predator[1]," (prey items=",length(bigpreyeaten),")",sep=""))
            
            lines(prop,BdistForPlotting,col="Black")
            firsttime = FALSE
          }else{
            lines(prop,BdistForPlotting,col="Black")
          }
        }
      }
      mtext("Contribution to diet",side=1,outer=TRUE,line=1.05)  
      mtext("Probability",side=2,outer=TRUE,line=1.05)
      dev.off() 
    }
  }
} #iloopends

npredlist  <- unique(x$Predator)


#open file

#write column names
columnnames<-c("PredatorName","PreyName","lower95","upper95", "mode", "simpleaverages", "bootaverage", "bootmode","normalizedmode")
#write to file here

preynames=names(x)[3:length(names(x))]

writeDF=columnnames
for(i in FirstPredToDo:LastPredToDo)  {
  
  mlesum = 0
  for(j in 1:nprey){
    mlesum = mlesum + mle[j,i]
  }
  
  for(j in 1:nprey){
    #newrow = c(toString(npredlist[i]),preynames[j],lower95[j,i],upper95[j,i],(mle[j,i]/mlesum),raw_wt_avg[j,i],bootaverage[j,i],bootmode[j,i])
    newrow = c(toString(npredlist[i]),preynames[j],lower95[j,i],upper95[j,i],mle[j,i],raw_wt_avg[j,i],bootaverage[j,i],bootmode[j,i],(mle[j,i]/mlesum))  
    writeDF=rbind(writeDF,newrow)
  }
}

write.csv(writeDF, "data/diet/Dirichlet/Dirichlet_Results.csv", row.names=FALSE)

