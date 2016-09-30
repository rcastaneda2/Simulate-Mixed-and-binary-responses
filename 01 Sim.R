

#Library ********************
library(mirt)
library(psych)
library(metafor)
library(mvtnorm)
library(logitnorm)
library(compiler)

load(".RData")

#Sources ********************
source("02 GenItems.R")
#source("Calculate Q3.R")

#p is items
#n is persons


#Moderators for a cluster of items
mod10.cluster<-c(rep(1,3),rep(0,42))
mod15.cluster<-c(rep(1,3),rep(0,102))
mod20.cluster<-c(rep(1,4),rep(0,186))
mod25.cluster<-c(rep(1,5),rep(0,295))
mod30.cluster<-c(rep(1,6),rep(0,429))
mod50.cluster<-c(rep(1,10),rep(0,1215))





#rho is the relationship between two thetas
sim<-function(reps,obs,items,ldmod,mix){#,rho,ldmod,complex,LD){
  p=items
  k= items*(items-1)/2
  v= rep(1/(obs-3),k)
  
  mod.simple = ldmod
  
  #If moderator specified then use that
  #if not then generate random moderator
  #Create matrix to capture results
  cell1out = matrix(nrow=reps, ncol=10)
  #Create list to capture parameter rmses 
  #rmse.res<-list(rmse="")
  #rmse.res<-list(rmse.res)[rep(1,reps)]
  if(mix=="binary"){
    rmse.res<-list(matrix(0,ncol=2,nrow=p),matrix(0,ncol=2,nrow=p))
  }else if(mix=="mix"){
    rmse.res<-list(matrix(0,ncol=4,nrow=p),matrix(0,ncol=4,nrow=p))
  }else if(mix=="graded"){
    rmse.res<-list(matrix(0,ncol=4,nrow=p),matrix(0,ncol=4,nrow=p))
  }else(print("Warning need binary, graded or TRUE"))
  
  i=1
  start=1
  while (i<reps){
    tryCatch( {
      print(c("restarted on ",i))
      
      start=i 
      for(i in start:reps){
          print(i)
          temp=genitems(obs,items,mix=mix)
          q3= .5*log((1+temp$q3)/(1-temp$q3))
          
          #RE simple
          out1=try(rma(q3~mod.simple,v,method='ML'))
            
          
          #CAPTURE DATA RE Simple #FULL
          cell1out[i,1] = out1$QE
          cell1out[i,2] = out1$QEp
          cell1out[i,3] = out1$QM
          cell1out[i,4] = out1$QMp #########################
          cell1out[i,5] = out1$b[1]
          cell1out[i,6] = out1$pval[1]
          cell1out[i,7] = out1$b[2]  ##################
          cell1out[i,8] = out1$pval[2]
          cell1out[i,9] = out1$tau2
          cell1out[i,10]= out1$R2
          
          rmse.res[[i]]<-temp$rmse
          
        }
        res=list(results=cell1out,rmse=rmse.res)
  return(res)
    },error=function(e){})
}

}



################ Real Simulation ################
################### CLUSTER #####################
#debug(sim)
#cell1<-sim(1000,1000,10,mod10.cluster,mix='mix')
#mean(cell1$results[,4]<.05)

#cell2<-sim(1000,1000,10,mod10.cluster,mix='binary')
#mean(cell2$results[,4]<.05)

#cell3<-sim(1000,1000,10,mod10.cluster,mix='graded')
#mean(cell3$results[,4]<.05)


cell13.1<-sim(1000,1000,10,mod10.cluster,mix='mix')
cell13.2<-sim(1000,1000,10,mod10.cluster,mix='mix')
cell13.3<-sim(1000,1000,10,mod10.cluster,mix='mix')
cell13.4<-sim(1000,1000,10,mod10.cluster,mix='mix')
cell13.5<-sim(1000,1000,10,mod10.cluster,mix='mix')
cell13.6<-sim(1000,1000,10,mod10.cluster,mix='mix'); save.image()


cell14.1<-sim(1000,1000,30,mod30.cluster,mix='mix')
cell14.2<-sim(1000,1000,30,mod30.cluster,mix='mix')
cell14.3<-sim(1000,1000,30,mod30.cluster,mix='mix')
cell14.4<-sim(1000,1000,30,mod30.cluster,mix='mix')
cell14.5<-sim(1000,1000,30,mod30.cluster,mix='mix')
cell14.6<-sim(1000,1000,30,mod30.cluster,mix='mix'); save.image()


cell15.1<-sim(1000,1000,50,mod50.cluster,mix='mix')
cell15.2<-sim(1000,1000,50,mod50.cluster,mix='mix')
cell15.3<-sim(1000,1000,50,mod50.cluster,mix='mix')
cell15.4<-sim(1000,1000,50,mod50.cluster,mix='mix')
cell15.5<-sim(1000,1000,50,mod50.cluster,mix='mix')
cell15.6<-sim(1000,1000,50,mod50.cluster,mix='mix'); save.image()



cell16.1<-sim(1000,4000,10,mod10.cluster,mix='mix')
cell16.2<-sim(1000,4000,10,mod10.cluster,mix='mix')
cell16.3<-sim(1000,4000,10,mod10.cluster,mix='mix')
cell16.4<-sim(1000,4000,10,mod10.cluster,mix='mix')
cell16.5<-sim(1000,4000,10,mod10.cluster,mix='mix')
cell16.6<-sim(1000,4000,10,mod10.cluster,mix='mix'); save.image()


cell17.1<-sim(1000,4000,30,mod30.cluster,mix='mix')
cell17.2<-sim(1000,4000,30,mod30.cluster,mix='mix')
cell17.3<-sim(1000,4000,30,mod30.cluster,mix='mix')
cell17.4<-sim(1000,4000,30,mod30.cluster,mix='mix')
cell17.5<-sim(1000,4000,30,mod30.cluster,mix='mix')
cell17.6<-sim(1000,4000,30,mod30.cluster,mix='mix'); save.image()


cell18.1<-sim(1000,4000,50,mod50.cluster,mix='mix')
cell18.2<-sim(1000,4000,50,mod50.cluster,mix='mix')
cell18.3<-sim(1000,4000,50,mod50.cluster,mix='mix')
cell18.4<-sim(1000,4000,50,mod50.cluster,mix='mix')
cell18.5<-sim(1000,4000,50,mod50.cluster,mix='mix')
cell18.6<-sim(1000,4000,50,mod50.cluster,mix='mix'); save.image()





















