

#data generating functions
#generate<-function(n,p,a,b,g,theta){
#  responses = matrix(0, nrow=n, ncol=p) 
#  for(i in 1:n) {
#    for(j in 1:p) {
#      prob =g[j]+(1-g[j])*((exp(a[j]*(theta[i]-b[j])))/(1+exp(a[j]*(theta[i]-b[j]))))
#      responses[i,j] = (runif(1) < prob)
#    }
#  }
#  return(responses)
#}

genitems<-function(n,p,mix='TRUE'){ #LD none vs cluster
  
  if (mix=='mix'){
    
            p1=floor(p*.3)
            p2=ceiling(p*.7)
            
            a=c(rnorm(p,1.7,.3))
            #First threshold
            b1=rnorm(p,-1.5,.5) #generate first b
            #second threshold
            b2=b1+rnorm(p,1,.2) #Generate second intercepts
            #Third threshold
            b3=b2+rnorm(p,1,.2)
            #Combine items
            b=cbind(b1,b2,b3)
            
            #Replace b for binary items
            binaryb<-rnorm(p2,0,1.5)
            b[(p1+1):p,1]<-binaryb
            #convert from b to d
            d<-(-1*(b/a))
            #replace thresholds on last n items to be binary
            d[(p1+1):p,2:3]<-NA
            items<-rep('dich',p)
            items[1:(p1)]<-'graded'
            
     
  }else if (mix=='binary'){
        #Only first two have LD
        
        a=c(rnorm(p,1.7,.3))
        #First threshold
        b=rnorm(p,0,1.5) #generate first b
        #second threshold
        #b2=b1+rnorm(n,1,.2) #Generate second intercepts
        #Third threshold
        #b3=b2+rnorm(n,1,.2)
        #Combine items
        #b=cbind(b1,b2,b3)
        #convert from b to d
        d<-(-1*(b/a))
        #replace thresholds on last n items to be binary
        #done[p1:p2,2:3]<-NA
        items<-rep('dich',p)
        
    }else if (mix=='graded'){
      #Only first two have LD
      
      a=c(rnorm(p,1.7,.3))
      #First threshold
      b1=rnorm(p,-1.5,.5) #generate first b
      #second threshold
      b2=b1+rnorm(p,1,.2) #Generate second intercepts
      #Third threshold
      b3=b2+rnorm(p,1,.2)
      #Combine items
      b=cbind(b1,b2,b3)
      #convert from b to d
      d<-(-1*(b/a))
      #replace thresholds on last n items to be binary
      #done[p1:p2,2:3]<-NA
      items<-rep('graded',p)
      
    }else( print ("ERROR!!"))
  
  #p1 = graded items
  #p2 = binary items
  #sample n number of times from normal distribution 

  
      theta<-matrix(rnorm(n))
      dat<-simdata(a,d,n,itemtype=items,Theta=theta)

      #Estimate parameters
      res=mirt(dat,1,verbose=FALSE)#,optimizer='L-BFGS-B',verbose=FALSE)
      temp3=residuals(res,type="Q3",verbose=FALSE)
      
      q3=matrix(temp3,ncol=p,nrow=p)
      
      #Parameters
      true.pars=cbind(a,d)
      
      est.pars=coef(res,simplify=TRUE)$items
      
      #Section on how close parameters are to estimated parameters. 
      
      
      if(mix=="mix"){
        colnames(true.pars)<-c('a1','d1','d2','d3')
        d3=rowSums(est.pars[,4:5],na.rm=TRUE)
        est.pars<-est.pars[,-c(4:7)]
        est.pars
        est.pars<-cbind(est.pars,d3)
        rmse<-(true.pars-est.pars)^2
        rmse<-round(rmse,4)
        
      }else if(mix=="binary"){
        colnames(true.pars)<-c('a1','d')
        est.pars<-est.pars[,-c(3:4)]  
        est.pars
        rmse<-(true.pars-est.pars)^2
        rmse<-round(rmse,4)
      }else if(mix=="graded"){
        colnames(true.pars)<-c('a1','d1','d2','d3')
        rmse<-(true.pars-est.pars)^2
        rmse<-round(rmse,4)
        
      }else(print("Error! Need TRUE, binary or graded"))
      
      
      #Set up matrix of Q3 to use in meta-analysis
      newq3=matrix(c(q3[upper.tri(q3,diag=FALSE)]),nrow=1)
      newq3=as.vector(newq3)
      
      

  final=list(dat=dat,q3=newq3,true.pars=true.pars,est.pars=est.pars,rmse=rmse)
  return(final)
  
}

#debug(genitems)
#debug(genitems)

#gen1<-genitems(4000,12,.2,LD="first")

undebug(genitems)
debug(genitems)
temp1<-genitems(1000,12,mix="binary")
#temp2<-genitems(1000,10,mix="binary")
#temp3<-genitems(1000,10,mix="graded")
