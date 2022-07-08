library(VGAM)   #fits Dirichlet Dist.
library(crayon) #color output to console
library(dplyr)
library(ggplot2)
library(viridis)
# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

### Update data, run Dirichlet ------------------------------------------------
run_Dirichlet <- function(data, name) {
  # Get total stomach weights for each predator and proportional weight for prey hake
  aged_wt <- data[, c("Predator_ID", "year", "predator_age", "prey_name", "prey_age", "prey_wt")] %>%
    group_by(Predator_ID) %>%
    mutate(stomach_wt = sum(prey_wt, na.rm = TRUE)) %>%
    filter(prey_name == "Pacific Hake") %>%
    mutate(hake_prey_prop = prey_wt / stomach_wt) %>%
    ungroup() %>%
    distinct()  # Remove duplicate rows - same pred ID, multiple hake prey
  
  aged_wt$new_ID <- c(1:nrow(aged_wt))
  
  # Subset dataset for only those hake that had no prey hake and assign proportion of diet as entirely "other"
  no_hake_dirichlet <- data %>%
    select(Predator_ID, year, predator_age, prey_name, prey_age, prey_wt) %>%
    group_by(Predator_ID) %>%
    filter(prey_name != "Pacific Hake") %>%
    ungroup() %>%
    select(Predator_ID, predator_age) %>%  
    mutate(predator_age = as.numeric(predator_age)) %>%
    distinct()
  no_hake_dirichlet$prey_age <- rep("other", nrow(no_hake_dirichlet))
  no_hake_dirichlet$hake_prey_prop <- rep(1, nrow(no_hake_dirichlet))
  
  # Select columns fron hake cannibalism dataset
  hake_dirichlet <- aged_wt %>%
    select(Predator_ID, predator_age, prey_age, hake_prey_prop)
  
  # Combine together and turn into wide format for Dirichlet script
  for_dirichlet <- rbind(hake_dirichlet, no_hake_dirichlet) %>%
    mutate(new_ID = c(1:length(Predator_ID))) %>%
    select(new_ID, predator_age, prey_age, hake_prey_prop) %>%
    pivot_wider(id_cols = c(new_ID, predator_age),  # rows to stay the same
                names_from = prey_age,  # columns to convert to wide
                values_from = hake_prey_prop,  # values to fill in
                values_fill = 0) %>%  # what to fill for missing values 
    arrange(predator_age)
  
  # Re-order and rename data for final dataset
  for_dirichlet <- for_dirichlet[, -1]
  colnames(for_dirichlet)[1] <- c("Predator")
  unique(as.character(for_dirichlet$Predator))
  for_dirichlet$Predator <- recode(for_dirichlet$Predator, 
                                   "1" = "pred_a1", "2" = "pred_a2",
                                   "3" = "pred_a3", "4" = "pred_a4",
                                   "5" = "pred_a5", "6" = "pred_a6",
                                   "7" = "pred_a7", "8" = "pred_a8",
                                   "9" = "pred_a9", "10" = "pred_a10",
                                   "11" = "pred_a11", "12" = "pred_a12",
                                   "13" = "pred_a13", "14" = "pred_a14",
                                   "15" = "pred_a15")
  
  # Add "wtg" column from Dirichlet example dataset
  for_dirichlet <- cbind(for_dirichlet[, 1], 
                         wtg = rep(1, length(for_dirichlet[, 1])),
                         for_dirichlet[, 2:(ncol(for_dirichlet) - 1)])
  
  # Update "other" column to include remainder from hake on hake stomachs
  for_dirichlet$other <- (1 - rowSums(for_dirichlet[, -c(1, 2), drop = FALSE]))
  
 
  ### Run Dirichlet script ------------------------------------------------------
  path <- "data/diet/Dirichlet/plots/"  # path for all plots
  x <- for_dirichlet
  
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
  
  ### Save output and plot results --------------------------------------------
  post_dirichlet <- data.frame(writeDF)
  
  # Fix issues from export from the Dirichlet script
  colnames(post_dirichlet) <- c("predator", "prey", "lower95", "upper95", "mode", 
                                "simple_average", "boot_average", "boot_mode",
                                "normalized_mode")
  post_dirichlet <- post_dirichlet[-1, ]
  
  post_dirichlet[, 3:9] <- sapply(post_dirichlet[, 3:9], as.numeric)  # make numeric
  
  post_dirichlet$predator <- factor(post_dirichlet$predator, 
                                    levels = c("pred_a1", "pred_a2", "pred_a3", 
                                               "pred_a4", "pred_a5", "pred_a6",
                                               "pred_a7", "pred_a8", "pred_a9", 
                                               "pred_a10", "pred_a11", "pred_a12",
                                               "pred_a13", "pred_a14", "pred_a15"))
  
  # Plot different statistics from the analysis
  post_dirichlet_long <- melt(post_dirichlet, id.vars = c("predator", "prey", "lower95", "upper95"))
  
  dirichlet_results <- ggplot(post_dirichlet_long, aes(x = prey)) +
    geom_errorbar(aes(ymin = lower95, ymax = upper95), 
                  width = .3, position = position_dodge(.9)) +
    geom_point(aes(x = prey, y = value, color = variable, shape = variable), size = 7, alpha = 0.5) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    xlab("prey item") + ylab("stomach proportion") +
    ggtitle(name) +
    facet_wrap(~ predator)
  
  return(list(post_dirichlet, dirichlet_results))
}

### Results -------------------------------------------------------------------
# Run for all years
# Read in full aged dataset
aged_dataset <- read.csv("data/diet/CCTD_FEAT_combined.csv")
# Replace NAs with age 1
aged_dataset$prey_age[is.na(aged_dataset$prey_age) & aged_dataset$prey_name == "Pacific Hake"] <- 1

all_years <- run_Dirichlet(aged_dataset, "All years")

all_years_df <- all_years[[1]]

all_years_plot <- all_years[[2]]
all_years_plot
ggsave(filename = "plots/diet/Dirichlet/Dirichlet_all_years.png", all_years_plot, 
       bg = "transparent", width=300, height=200, units="mm", dpi=300)

# Run for only 1998
y90s <- c(1991, 1995, 1997, 1998, 1999)
df_90s <- aged_dataset %>% filter(year %in% y90s)
dirichlet_90s <- run_Dirichlet(df_90s, "1991-1999")
dirichlet_90s[[2]]

# Run with recent data
recent <- c(2005, 2007, 2011, 2015, 2019)
df_recent <- aged_dataset %>% filter(year %in% recent)
dirichlet_recent <- run_Dirichlet(df_recent, "2005-2019")
dirichlet_recent[[2]]


### Re-organize dataframe for CEATTLE -----------------------------------------
to_ceattle <- function(df) {
  dirichlet_ceattle <- df %>% filter(prey != "other")
  
  dirichlet_ceattle <- data.frame(cbind(Pred = rep(1, nrow(dirichlet_ceattle)),
                                        Prey = rep(1, nrow(dirichlet_ceattle)),
                                        Pred_sex = rep(0, nrow(dirichlet_ceattle)),
                                        Prey_sex = rep(0, nrow(dirichlet_ceattle)),
                                        Pred_age = dirichlet_ceattle$predator,
                                        Prey_age = dirichlet_ceattle$prey,
                                        Year = rep(0, nrow(dirichlet_ceattle)),
                                        Sample_size = rep(10, nrow(dirichlet_ceattle)),
                                        Stomach_proportion_by_weight = dirichlet_ceattle$boot_average))
}

ceattle_all <- to_ceattle(all_years_df)
write.csv(ceattle_all, "data/diet/Dirichlet/Dirichlet_all_years.csv", row.names = FALSE)

ceattle_90s <- to_ceattle(dirichlet_90s[[1]])
write.csv(ceattle_90s, "data/diet/Dirichlet/Dirichlet_90s.csv", row.names = FALSE)

ceattle_recent <- to_ceattle(dirichlet_recent[[1]])
write.csv(ceattle_recent, "data/diet/Dirichlet/Dirichlet_recent.csv", row.names = FALSE)
