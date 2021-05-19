library(RANN)


CASOS<-function(X,X.weight=NULL,weight.scheme=NULL,np=3,ASmethod="LOF",k=5,h=1,standardise=TRUE,do.all.dim=FALSE){
	# X = data matrix;
	#     all columns should be magnitude columns if used with astro data
      # X.weight = matrix with weights per variable
      #            will be used to compute the weights to apply to each subspace
      #            if NULL then all weights are equal and set to 1
      # weight.scheme = how to use X.weight to weight the anomaly score (one of "ratio1","ratio2","ratio0-5","ratioAVG", "exp1","exp2","exp0-5" cf. further below in code for more info)
	# np = maximum dimensionality of projections;
	#      all projections of all dimensinalities<=np will be computed if do.al.dim==T, otherwise only projections of dimensionality exactly np will be computed
	# ASmethod = method to be used to compute anomaly scores;
	#            for the moment are supported:
	#               "LOF" (Local Outlier Factor; Breuning et al. ACM SIGMOD Conf. Proc.; 2000)
	#               "LDF" (Outlier Detection with Kernel Density Functions; Latecki, Lazarevic and Pokrajac, LNAI 4571; 2007)
	#               "LDF.ecld" (a version of LDF using only Euclidean distance rather than Euclidean to find the kNN and then Mahalanobis for score calculation)
	# k = number of nearest neighbours (for ASmethod="LOF", "LDF" and "LDF.ecld")
	# h = smoothing parameter (for ASmethod="LDF" and "LDF.ecld")
	# standardise = logical;
	#               if TRUE all columns are standardised to zero mean and unit variance
	# required package: dprep (if ASmethod="LOF")
      # do.all.dim = logical;
      #         if TRUE then AS in all subspaces with dim<=np will be computed, else only subspaces of dimesnionality np will be computed
	
	
	nc<-ncol(X)
	nr<-nrow(X)
	if(nc<np){np<-nc; print("ncol(X) < np. np has been reset to ncol(X).")}
	if(standardise==TRUE){
		for(i in 1:nc){
			if(sum(is.na(X[,i]))<nr){X[,i]<-(X[,i]-mean(X[!is.na(X[,i]),i]))/sd(X[!is.na(X[,i]),i])}
		}
	}
		
	## FINDPROJECTIONS AND COMPUTE ANOMALY SCORES
		
		# projList will consist of np lists, each consisting of combinatorial(n=nc,k=k) lists
		#          where k refers to the current list of projList
		# projList contains the column indices corresponding to the different projections that will be used
		projList<-vector("list")
		
		NP<-as.integer(0)
		  # NP = total number of different projections to compute
            LOF.sc<-vector("list")
		  # LOF.sc will become the output matrix containing the various LOF scores
		
		if(do.all.dim==T){mindim<-1}
		if(do.all.dim==F){mindim<-np}
		
		for(i in mindim:np){

			NP<-as.integer(NP+round(as.integer(combinatorial(n=nc,k=i))))
			projList[[i]]<-vector("list")
                  LOF.sc[[i]]<-matrix(numeric(0),nrow=nrow(X))
			check<-as.integer(1)

			while(check<=as.integer(round(combinatorial(nc,i)))){
				
				# find a subspace not used in a previous iteration
				  projList[[i]][[check]]<-sort(sample(i,x=1:nc,replace=FALSE))
				  if(length(unique(projList[[i]]))==check){
				  	
				  	# find all data objects that have observations in the selected subspace
				  	  Xtemp<-X[,projList[[i]][[check]]]
                                if(!is.null(X.weight) & i>1){weights<-sqrt(apply(X=(X.weight[,projList[[i]][[check]]])^2,MARGIN=1,FUN=sum))}
                                if(!is.null(X.weight) & i==1){weights<-sqrt(apply(X=(as.matrix(X.weight[,projList[[i]][[check]]],ncol=1))^2,MARGIN=1,FUN=sum))}
				  	  if(length(projList[[i]][[check]])>1){
				  	  	#idx<-vector("numeric",nrow(X))
				  	  	#for(j in 1:ncol(Xtemp)){idx<-idx+Xtemp[,j]}
				  	  	idx<-vector("integer",nrow(X))
                                    for(j in 1:ncol(Xtemp)){idx[is.na(Xtemp[,j])]<-NA}
				  	  }
				  	  if(length(projList[[i]][[check]])==1){
				  	  	#idx<-vector("numeric",nrow(X))
				  	  	#idx<-Xtemp
				  	  	idx<-vector("integer",nrow(X))
				  	  	idx[is.na(Xtemp)]<-NA
				  	  }

				  	# compute the anomaly scores
				  	
				  	  if(ASmethod=="LOF"){
				  	  	#library(dprep) # needed if the lofactor subroutine is used
				  	  	#not used now, since we have coded a more computationally efficient version - lofactor.opt
				  	    LOF.sc[[i]]<-cbind(LOF.sc[[i]],rep(numeric(0),length=nrow(X)))
				  	    if(sum(!is.na(idx))>k & length(projList[[i]][[check]])>1){LOF.sc[[i]][!is.na(idx),check]<-lofactor.opt(data=Xtemp[!is.na(idx),],k=k)}
				  	    if(sum(!is.na(idx))>k & length(projList[[i]][[check]])==1){LOF.sc[[i]][!is.na(idx),check]<-lofactor.opt(data=Xtemp[!is.na(idx)],k=k)}
				  	    if(sum(!is.na(idx))<=k & sum(!is.na(idx))>2 & length(projList[[i]][[check]])>1){LOF.sc[[i]][!is.na(idx),check]<-lofactor.opt(data=Xtemp[!is.na(idx),],k=(sum(!is.na(idx))-1))}
				  	    if(sum(!is.na(idx))<=k & sum(!is.na(idx))>2 & length(projList[[i]][[check]])==1){LOF.sc[[i]][!is.na(idx),check]<-lofactor.opt(data=Xtemp[!is.na(idx)],k=(sum(!is.na(idx))-1))}
				  	    if(sum(!is.na(idx))<=2){LOF.sc[[i]][,check]<-NA}
				  	    rm(Xtemp)
				  	  }
				  	  if(ASmethod=="LDF"){
				  	    LOF.sc[[i]]<-cbind(LOF.sc[[i]],rep(numeric(0),length=nrow(X)))
				  	    if(sum(!is.na(idx))>k & length(projList[[i]][[check]])>1){LOF.sc[[i]][!is.na(idx),check]<-ldfactor(data=Xtemp[!is.na(idx),],k=k,h=h)}
				  	    if(sum(!is.na(idx))>k & length(projList[[i]][[check]])==1){LOF.sc[[i]][!is.na(idx),check]<-ldfactor(data=Xtemp[!is.na(idx)],k=k,h=h)}
				  	    if(sum(!is.na(idx))<=k & sum(!is.na(idx))>2 & length(projList[[i]][[check]])>1){LOF.sc[[i]][!is.na(idx),check]<-ldfactor(data=Xtemp[!is.na(idx),],k=(sum(!is.na(idx))-1),h=h)}
				  	    if(sum(!is.na(idx))<=k & sum(!is.na(idx))>2 & length(projList[[i]][[check]])==1){LOF.sc[[i]][!is.na(idx),check]<-ldfactor(data=Xtemp[!is.na(idx)],k=(sum(!is.na(idx))-1),h=h)}
				  	    if(sum(!is.na(idx))<=2){LOF.sc[[i]][,check]<-NA}
				  	    rm(Xtemp)
				  	  }
				  	  if(ASmethod=="LDF.ecld"){
				  	    LOF.sc[[i]]<-cbind(LOF.sc[[i]],rep(numeric(0),length=nrow(X)))
				  	    if(sum(!is.na(idx))>k & length(projList[[i]][[check]])>1){LOF.sc[[i]][!is.na(idx),check]<-ldfactor.ecld(data=Xtemp[!is.na(idx),],k=k,h=h)}
				  	    if(sum(!is.na(idx))>k & length(projList[[i]][[check]])==1){LOF.sc[[i]][!is.na(idx),check]<-ldfactor.ecld(data=Xtemp[!is.na(idx)],k=k,h=h)}
				  	    if(sum(!is.na(idx))<=k & sum(!is.na(idx))>2 & length(projList[[i]][[check]])>1){LOF.sc[[i]][!is.na(idx),check]<-ldfactor.ecld(data=Xtemp[!is.na(idx),],k=(sum(!is.na(idx))-1),h=h)}
				  	    if(sum(!is.na(idx))<=k & sum(!is.na(idx))>2 & length(projList[[i]][[check]])==1){LOF.sc[[i]][!is.na(idx),check]<-ldfactor.ecld(data=Xtemp[!is.na(idx)],k=(sum(!is.na(idx))-1),h=h)}
				  	    if(sum(!is.na(idx))<=2){LOF.sc[[i]][,check]<-NA}
				  	    rm(Xtemp)
				  	  }
				  	            ## different weighting schemes: either w=a/(a+err) with a=... ("ratio1,2,etc") or w=exp(-a*(err^2)/2) with a=... ("exp1,2,etc")
                                if(!is.null(X.weight)){
                                	if(weight.scheme=="ratio1"){LOF.sc[[i]][!is.na(idx),check]<-LOF.sc[[i]][!is.na(idx),check]/(1+weights[!is.na(idx)])}
                                	if(weight.scheme=="ratio2"){LOF.sc[[i]][!is.na(idx),check]<-2*LOF.sc[[i]][!is.na(idx),check]/(2+weights[!is.na(idx)])}
                                	if(weight.scheme=="ratio0-5"){LOF.sc[[i]][!is.na(idx),check]<-0.5*LOF.sc[[i]][!is.na(idx),check]/(0.5+weights[!is.na(idx)])}
                                	if(weight.scheme=="ratioAVG"){m.temp<-mean(weights[!is.na(idx)]); LOF.sc[[i]][!is.na(idx),check]<-m.temp*LOF.sc[[i]][!is.na(idx),check]/(m.temp+weights[!is.na(idx)])}
                                	if(weight.scheme=="exp1"){LOF.sc[[i]][!is.na(idx),check]<-LOF.sc[[i]][!is.na(idx),check]*exp(-(weights[!is.na(idx)]^2)/2)}
                                	if(weight.scheme=="exp2"){LOF.sc[[i]][!is.na(idx),check]<-LOF.sc[[i]][!is.na(idx),check]*exp(-2*(weights[!is.na(idx)]^2)/2)}
                                	if(weight.scheme=="exp0-5"){LOF.sc[[i]][!is.na(idx),check]<-LOF.sc[[i]][!is.na(idx),check]*exp(-(0.5)*(weights[!is.na(idx)]^2)/2)}
                                }

                            check<-as.integer(round(check+as.integer(1)))
				  }	  
			}
		}	
	res<-vector("list",4)
	names(res)<-c("LOF.sc","NP","projList","ASmethod")
	res$LOF.sc<-LOF.sc
	res$NP<-NP
	res$projList<-projList
	res$ASmethod<-ASmethod
	rm(LOF.sc,NP,projList)
	res
}


# auxiliary sub-routines


lofactor.opt<-function(data,k,standardise=F){
  ## data = data matrix with rows = objects, columns = variables
  ## k = number of nearest neighbours to consider
  ## N.B.: nn2 cannot handle ties for nearest neighbour!!!

  #library(RANN) # required for the nn2 subroutine
  
  if(standardise==TRUE){
  	n<-nrow(data)
  	p<-ncol(data)
    for(i in 1:p){
      if(sum(is.na(data[,i]))<n){data[,i]<-(data[,i]-mean(data[!is.na(data[,i]),i]))/sd(data[!is.na(data[,i]),i])}
    }
  }

  # nearest neighbour indices, k-distances
    if(!is.null(nrow(data))){
      n<-nrow(data)
      p<-ncol(data)
      nn.list<-nn2(data=data,query=data,k=(k+1),eps=0)
      kdist<-nn.list$nn.dists[,k+1]
    }
    if(is.null(nrow(data))){
      n<-length(data)
      p<-1
      data.opt<-cbind(data,data)
      nn.list<-nn2(data=data.opt,query=data.opt,k=(k+1),eps=0)
      for(i in 2:(k+1)){
        nn.list$nn.dists[,i]<-nn.list$nn.dists[,i]/sqrt(2)
      }
      kdist<-nn.list$nn.dists[,k+1]
      rm(data.opt)
    }

  # local reachability density
    lrd<-numeric(0)
    for(i in 1:n){
      temp1<-kdist[nn.list$nn.idx[i,c(2:(k+1))]]
      temp2<-nn.list$nn.dist[i,c(2:(k+1))]
      rd<-temp1
      rd[temp1<temp2]<-temp2[temp1<temp2]
      lrd[i]<-1/(mean( rd ))
    }

  # LOF
    Xtemp<-matrix(lrd[nn.list$nn.idx[1:n,-1]],ncol=k,byrow=F)
    lofval<-as.vector(apply(X=Xtemp,MARGIN=1,FUN=mean))
    lofval<-lofval/lrd
    rm(Xtemp)

  lofval

}


Comb.Ext<-function(ASmat){
    n<-nrow(ASmat)
    AD.Ext<-vector("numeric",n)
	for(j in 1:n){
        lof.sc<-ASmat[j,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        AD.Ext[j]<-max(lof.sc)
    }
    rm(n)
    AD.Ext
}


Comb.Ave<-function(ASmat){
    n<-nrow(ASmat)
    AD.Ave<-vector("numeric",n)
	for(j in 1:n){
        lof.sc<-ASmat[j,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        AD.Ave[j]<-mean(lof.sc)
    }
    rm(n)
    AD.Ave    
}


Comb.Med<-function(ASmat,pQuant=0.5){
    n<-nrow(ASmat)
    AD.Med<-vector("numeric",n)
	for(j in 1:n){
        lof.sc<-ASmat[j,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        AD.Med[j]<-quantile.opt(lof.sc,probs=pQuant)
    }
    rm(n)
    AD.Med    
}


Comb.TopN<-function(ASmat,N){
    if(N>ncol(ASmat)){N<-ncol(ASmat); print("N was too large. Reduced to maximum -- TopN now equals Ave.")}
    n<-nrow(ASmat)
    AD.TopN<-vector("numeric",n)
	for(j in 1:n){
        lof.sc<-ASmat[j,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        Nalt<-N
	  nalt<-length(lof.sc)
        if(N>nalt){Nalt<-nalt}
        AD.TopN[j]<-mean((sort(lof.sc))[(nalt-Nalt+1):nalt])
        rm(Nalt,nalt)
    }
    rm(n)
    AD.TopN
}


Comb.OneBadVar<-function(ASmat,np){
    n<-nrow(ASmat)
    AD.OneBadVar<-vector("numeric",n)
    for(i in 1:n){
      s<-sum(!is.na(ASmat[i,]))
      if(s>=np){
        # determinging the number nc of observed variables for the given object
          for(j in np:s){
            if(s==factorial(j)/(factorial(np)*(factorial(j-np)))){nc<-j}
          }
        # if there is one badly measured variable, it will affect combinatorial((nc-1),(np-1)) ASs
          t<-combinatorial((nc-1),(np-1))
          AD.OneBadVar[i]<-sort(ASmat[i,],decreasing=T,na.last=T)[t+1]
      }
      if(s<np){AD.OneBadVar[i]<-NA}
    }
    rm(s,t,nc,n)
    AD.OneBadVar
}


Comb.TopTail<-function(ASmat,pTail){
    n<-nrow(ASmat)
    AD.TopTail<-vector("numeric",n)
	TopTail.val<-apply(X=ASmat,MARGIN=2,FUN=quantile.opt,probs=pTail)
	for(j in 1:n){
        lof.sc<-ASmat[j,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        TopTail.val.alt<-TopTail.val[!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        lof.sc<-lof.sc-TopTail.val.alt
        lof.sc[lof.sc<0]<-0
        AD.TopTail[j]<-sum(lof.sc)
        rm(lof.sc,TopTail.val.alt)
    }
    rm(TopTail.val,n)
    AD.TopTail
}


Comb.Mid<-function(ASmat,pMid){
    if(pMid[1]>pMid[2]){stop("Need pMid[1]<=pMid[2].")}
    n<-nrow(ASmat)
    AD.Mid<-vector("numeric",n)
	Mid.val<-apply(X=ASmat,MARGIN=2,FUN=quantile.opt,probs=pMid)
	for(j in 1:n){
        lof.sc<-ASmat[j,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        Mid.val.low<-Mid.val[1,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        Mid.val.high<-Mid.val[2,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        lof.sc[lof.sc>Mid.val.high]<-0
        lof.sc<-lof.sc-Mid.val.low
        lof.sc[lof.sc<0]<-0
        AD.Mid[j]<-sum(lof.sc)
        rm(lof.sc,Mid.val.low,Mid.val.high)
    }
    rm(Mid.val,n)
    AD.Mid
}


Comb.MidSoft<-function(ASmat,pMid){
    if(pMid[1]>pMid[2]){stop("Need pMid[1]<=pMid[2].")}
    n<-nrow(ASmat)
    AD.Mid<-vector("numeric",n)
	Mid.val<-apply(X=ASmat,MARGIN=2,FUN=quantile.opt,probs=pMid)
	for(j in 1:n){
        lof.sc<-ASmat[j,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        Mid.val.low<-Mid.val[1,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        Mid.val.high<-Mid.val[2,!is.na(ASmat[j,]) & abs(ASmat[j,])<Inf]
        #for(s in 1:length(lof.sc)){if(lof.sc[s]>Mid.val.high[s]){lof.sc[s]<-Mid.val.high[s]}}
        #for(s in 1:length(lof.sc)){if(lof.sc[s]>Mid.val.high[s]){lof.sc[s]<-Mid.val.high[s]-(lof.sc[s]-Mid.val.high[s])}}
        lof.sc[lof.sc>Mid.val.high]<-(Mid.val.high-(lof.sc-Mid.val.high))[lof.sc>Mid.val.high]
        lof.sc<-lof.sc-Mid.val.low
        lof.sc[lof.sc<0]<-0
        AD.Mid[j]<-sum(lof.sc)
        rm(lof.sc,Mid.val.low,Mid.val.high)
    }
    rm(Mid.val,n)
    AD.Mid
}


CASOS.combine<-function(ASmat,combfun,...){
  # combines AS
  # it is applicable to individual AS matrices
  # There are 8 different combination scheme to be specified as combfun (each can be specified in two ways):
  #    "Ext" / Comb.Ext: taking the most extreme ASs,
  #    "Ave" / Comb.Ave: averaging all ASs
  #    "Med" / Comb.Med: taking the median of all ASs
  #    "TopN" / Comb.TopN: averaging the top N ASs
  #    "TopTail" / Comb.TopTail: summing the ASs above a certain percentile minus that percentile per object
  #    "Mid" / Comb.Mid: summing the ASs between two percentiles minus the lower percentile per object
  #    "MidSoft" / Comb.MidSoft: as "Mid" above, but ASs above the higher percentile are set to that percentile
  #    "OneBadVar" / Comb.OneBadVar: in the case there is one badly measured variable, we will take as AS the Nth highest AS, where N is the number of affected subspaces
  # ASmat = an object of type AnomDet (i.e. output from the AnomDet routine; a matrix with AS)
  # N     = number of top ASs scores to average in the Top-N combination scheme
  # pTail = percentile above which ASs have to be to be taken into consideration for TopTail scheme
  # pMid = percentiles between which the ASs have to be taken into consideration for Mid & MidSoft schemes
  # combfun = cobination function; one of "Comb.Ext", "Comb.Ave", "Comb.TopN", "Comb.TopTail"
  # ... = any parameters required for combfun
  
    p<-ncol(ASmat)
    n<-nrow(ASmat)
    
    if(!is.function(combfun)){if(combfun=="Ext"){combfun<-Comb.Ext}}
    if(!is.function(combfun)){if(combfun=="Ave"){combfun<-Comb.Ave}}
    if(!is.function(combfun)){if(combfun=="Med"){combfun<-Comb.Med}}
    if(!is.function(combfun)){if(combfun=="TopN"){combfun<-Comb.TopN}}
    if(!is.function(combfun)){if(combfun=="TopTail"){combfun<-Comb.TopTail}}
    if(!is.function(combfun)){if(combfun=="Mid"){combfun<-Comb.Mid}}
    if(!is.function(combfun)){if(combfun=="MidSoft"){combfun<-Comb.MidSoft}}
    if(!is.function(combfun)){if(combfun=="OneBadVar"){combfun<-Comb.OneBadVar}}

    AD.comb<-combfun(ASmat,...)
    AD.comb    
}


ldfactor<-function(data,k,h,standardise=F){
  ## ldf score computation as in Latecki, Lazarevic and Pokrajac (2007), LNAI 4571, pp.61-75
  ## data = data matrix with rows = objects, columns = variables
  ## k = number of nearest neighbours to consider
  ## h = smoothing factor
  ## N.B.: nn2 cannot handle ties for nearest neighbour!!!

  #library(RANN) # required for the nn2 subroutine
  
  if(standardise==TRUE){
  	n<-nrow(data)
  	p<-ncol(data)
    for(i in 1:p){
      if(sum(is.na(data[,i]))<n){data[,i]<-(data[,i]-mean(data[!is.na(data[,i]),i]))/sd(data[!is.na(data[,i]),i])}
    }
  }

  # nearest neighbour indices, k-distances
    if(!is.null(nrow(data))){
      n<-nrow(data)
      p<-ncol(data)
      nn.list<-nn2(data=data,query=data,k=(k+1),eps=0)
      kdist<-nn.list$nn.dists[,k+1]
    }
    if(is.null(nrow(data))){
      n<-length(data)
      p<-1
      data.opt<-cbind(data,data)
      nn.list<-nn2(data=data.opt,query=data.opt,k=(k+1),eps=0)
      for(i in 2:(k+1)){
        nn.list$nn.dists[,i]<-nn.list$nn.dists[,i]/sqrt(2)
      }
      kdist<-nn.list$nn.dists[,k+1]
      rm(data.opt)
    }

  # local density estimate
    lde<-numeric(0)
    for(i in 1:n){
      rd<-numeric(0)
      dist.to.knn<-rep(0,k+1)
      det.sigma<-numeric(0)
      for(j in 2:(k+1)){
        data.temp<-data[nn.list$nn.idx[nn.list$nn.idx[i,j],],]
        Sigma<-var(data.temp)
        det.sigma[j-1]<-det(Sigma)
        mu<-apply(X=data.temp,MARGIN=2,FUN=mean)
        eigen.dec<-eigen(Sigma,symmetric=T,only.values=F)
        data.alt<-t( (sqrt(diag(1/eigen.dec$values)) %*% t(eigen.dec$vectors)) %*% (t(data.temp)-mu) )
        #Sigma1<-solve(Sigma)
        dist.temp<-rep(0,k+1)
        for(s in 1:(k+1)){
          dist.temp[s]<-sqrt( sum( (data.alt[1,]-data.alt[s,])^2 ) )
          #dist.temp[s]<-sqrt( ((data.temp[1,]-data.temp[s,])%*%Sigma1)%*%(data.temp[1,]-data.temp[s,]) )
        }
        dist.to.knn[j]<-max(dist.temp)
        rm(data.temp,dist.temp,data.alt)
      }
      data.temp<-data[nn.list$nn.idx[i,],]
      Sigma<-var(data.temp)
      mu<-apply(X=data.temp,MARGIN=2,FUN=mean)
      eigen.dec<-eigen(Sigma,symmetric=T,only.values=F)
      data.alt<-t( (sqrt(diag(1/eigen.dec$values)) %*% t(eigen.dec$vectors)) %*% (t(data.temp)-mu) )
      #Sigma1<-solve(Sigma)
      for(j in 1:k){
        temp<-dist.to.knn[j+1]
        if(p>1){rd[j]<-max( temp,sqrt( sum((data.alt[j+1,]-data.alt[1,])^2) ) )}
        if(p==1){rd[j]<-max( temp,sqrt( sum((data.alt[j+1]-data.alt[1])^2) ) )}
        #if(p>1){rd[j]<-max( temp,sqrt( ((data.temp[1,]-data.temp[j+1,])%*%Sigma1)%*%(data.temp[1,]-data.temp[j+1,]) ) )}
        #if(p==1){rd[j]<-max( temp,sqrt( ((data.temp[1]-data.temp[j+1])^2)*Sigma1 ) )}
      }
      lde[i]<-mean( (1/(sqrt(det.sigma)*(sqrt(2*pi)*h)^p)) * exp( -rd^2/(2*h^2) ) )
    }
    rm(temp,data.temp,data.alt,dist.to.knn,rd)
    
  # LDF
    Xtemp<-matrix(lde[nn.list$nn.idx[1:n,-1]],ncol=k,byrow=F)
    ldfval<-as.vector(apply(X=Xtemp,MARGIN=1,FUN=mean))
    ldfval<-ldfval/(lde+0.1*ldfval)
    #ldfval<-ldfval/lde
    rm(Xtemp)

  ldfval
}


ldfactor.ecld<-function(data,k,h,standardise=F){
  ## as ldfactor above, just with Euclidean distance used throughout
  ## data = data matrix with rows = objects, columns = variables
  ## k = number of nearest neighbours to consider
  ## h = smoothing factor
  ## N.B.: nn2 cannot handle ties for nearest neighbour!!!

  #library(RANN) # required for the nn2 subroutine
  
  if(standardise==TRUE){
  	n<-nrow(data)
  	p<-ncol(data)
    for(i in 1:p){
      if(sum(is.na(data[,i]))<n){data[,i]<-(data[,i]-mean(data[!is.na(data[,i]),i]))/sd(data[!is.na(data[,i]),i])}
    }
  }

  # nearest neighbour indices, k-distances
    if(!is.null(nrow(data))){
      n<-nrow(data)
      p<-ncol(data)
      nn.list<-nn2(data=data,query=data,k=(k+1),eps=0)
      kdist<-nn.list$nn.dists[,k+1]
    }
    if(is.null(nrow(data))){
      n<-length(data)
      p<-1
      data.opt<-cbind(data,data)
      nn.list<-nn2(data=data.opt,query=data.opt,k=(k+1),eps=0)
      for(i in 2:(k+1)){
        nn.list$nn.dists[,i]<-nn.list$nn.dists[,i]/sqrt(2)
      }
      kdist<-nn.list$nn.dists[,k+1]
      rm(data.opt)
    }

  # local density estimate
    lde<-numeric(0)
    for(i in 1:n){
      rd<-numeric(0)
      det.sigma<-numeric(0)  
      for(j in 1:k){
      	data.temp<-data[nn.list$nn.idx[nn.list$nn.idx[i,j+1],],]
      	det.sigma[j]<-det(var(data.temp))
      }
      temp1<-kdist[nn.list$nn.idx[i,c(2:(k+1))]]
      temp2<-nn.list$nn.dist[i,c(2:(k+1))]
      rd<-temp1
      rd[temp1<temp2]<-temp2[temp1<temp2]
      lde[i]<-mean( (1/((sqrt(det.sigma)*(sqrt(2*pi)*h)^p))) * exp( -rd^2/(2*h^2) ) )
    }
    #rm(temp)
     
  # LDF
    Xtemp<-matrix(lde[nn.list$nn.idx[1:n,-1]],ncol=k,byrow=F)
    ldfval<-as.vector(apply(X=Xtemp,MARGIN=1,FUN=mean))
    ldfval<-ldfval/(lde+0.1*ldfval)
    #ldfval<-ldfval/lde
    rm(Xtemp)

  ldfval
}


loci<-function(data,ksigma=3,nmin=20,alpha=0.5,standardise=F,verbose=F){
  ## computes LOCI score as in Papadimitriou, Kitagawa and Gibbons (2003), Proc. of ICDE'03, p.315-323
  ## data = data matrix with rows = objects, columns = variables
  ## ksigma = threshold factor for cutting off anomalies (i.e. how many stds)
  ## nmin = minimum number of NNs within a given radius
  ## alpha = by how much to extend neighbourhood? (alpha<1)
  ## method is harsh on memory requirements due to needing to store the entire distance matrix while it is running
  ## incredibly slow at the moment - not workable

  n<-nrow(data)
  flag<-integer(n)


  # pre-processing
    if(standardise==TRUE){
      n<-nrow(data)
      p<-ncol(data)
      for(i in 1:p){
        if(sum(is.na(data[,i]))<n){data[,i]<-(data[,i]-mean(data[!is.na(data[,i]),i]))/sd(data[!is.na(data[,i]),i])}
      }
    }

    D<-as.matrix(dist(data,method="euclidean"))

    if(verbose){print("Pre-processing done.")}

  # post-processing
    for(i in 1:n){
      #rmin<-sort(as.matrix(D)[i,])[nmin]
      dist.temp<-D[i,]
      r.temp<-sort(dist.temp)[-(1:nmin)]
      for(r in r.temp){
        if(verbose){print(paste("SAMPLE = ",i,", radius = ",r,sep=""))}
        n_alphar<-sum( dist.temp<alpha*r )
        idx.NN<-(1:n)[dist.temp<r]
        n_alphar_NNs<-integer(length(idx.NN))
        for(j in 1:length(idx.NN)){n_alphar_NNs[j]<-sum( as.matrix(D)[idx.NN[j],]<alpha*r )}
        nhat_alphar<-mean(n_alphar_NNs)
        MDEF<-(1-n_alphar/nhat_alphar)
        sigmaMDEF<-sqrt( sum( (n_alphar_NNs-rep(nhat_alphar,length=length(idx.NN)))^2 )/length(idx.NN) )/nhat_alphar
        #print(paste("..",n_alphar,length(idx.NN),nhat_alphar,MDEF,sigmaMDEF))
        if(MDEF>ksigma*sigmaMDEF){flag[i]<-1}
      }
    }
  flag
}


loci.opt<-function(data,ksigma=3,nmin=20,alpha=0.5,nrad=10,standardise=F,verbose=F){
  ## computes LOCI score as in Papadimitriou, Kitagawa and Gibbons (2003), Proc. of ICDE'03, p.315-323
  ## data = data matrix with rows = objects, columns = variables
  ## ksigma = threshold factor for cutting off anomalies (i.e. how many stds)
  ## nmin = minimum number of NNs within a given radius
  ## alpha = by how much to extend neighbourhood? (alpha<1)
  ## method is harsh on memory requirements due to needing to store the entire distance matrix while it is running
  ## incredibly slow at the moment - not workable

  n<-nrow(data)
  flag<-integer(n)

  minmax.finder<-function(x,nmin){
    sort(x)[c(nmin,length(x))]
  }

  # pre-processing
    if(standardise==TRUE){
      n<-nrow(data)
      p<-ncol(data)
      for(i in 1:p){
        if(sum(is.na(data[,i]))<n){data[,i]<-(data[,i]-mean(data[!is.na(data[,i]),i]))/sd(data[!is.na(data[,i]),i])}
      }
    }
  
    D<-as.matrix(dist(data,method="euclidean"))
    dist.temp<-apply(X=D,FUN=minmax.finder,MARGIN=1,nmin=nmin) # any distance can be used; here we use Euclidean for now
    dist.min<-min(dist.temp[1,])
    dist.max<-max(dist.temp[2,])
    rm(dist.temp,minmax.finder)

    if(verbose){print("Pre-processing done.")}

  # post-processing
    for(i in 1:n){
      #rmin<-sort(as.matrix(D)[i,])[nmin]
      dist.temp<-D[i,]
      for(r in seq(from=dist.min,to=dist.max,length=nrad)){
        if(verbose){print(paste("SAMPLE = ",i,", radius = ",r,sep=""))}
        n_alphar<-sum( dist.temp<alpha*r )
        idx.NN<-(1:n)[dist.temp<r]
        n_alphar_NNs<-integer(length(idx.NN))
        for(j in 1:length(idx.NN)){n_alphar_NNs[j]<-sum( as.matrix(D)[idx.NN[j],]<alpha*r )}
        nhat_alphar<-mean(n_alphar_NNs)
        MDEF<-(1-n_alphar/nhat_alphar)
        sigmaMDEF<-sqrt( sum( (n_alphar_NNs-rep(nhat_alphar,length=length(idx.NN)))^2 )/length(idx.NN) )/nhat_alphar
        #print(paste("..",n_alphar,length(idx.NN),nhat_alphar,MDEF,sigmaMDEF))
        if(MDEF>ksigma*sigmaMDEF){flag[i]<-1}
      }
    }
  flag
}


abod.outlier<-function(dat,standardise=F,verbose=F){
  ## method described in Kriegel et al (2008)

  nc<-ncol(dat)
  nr<-nrow(dat)

  if(standardise==TRUE){
    for(i in 1:nc){
      if(sum(is.na(dat[,i]))<nr){dat[,i]<-(dat[,i]-mean(dat[!is.na(dat[,i]),i]))/sd(dat[!is.na(dat[,i]),i])}
    }
  }

  AS<-rep(NA,nr)
  for(i in 1:nr){
    angles.vect<-numeric(0)
    for(j1 in (1:nr)[-i]){
      for(j2 in (1:nr)[-c(i,j1)]){
        if(verbose){print(paste("SAMPLE = ",i,", j1 = ",j1,", j2 = ",j2,sep=""))}
        angles.vect<-c(  angles.vect,sum( (dat[j1,]-dat[i,])*(dat[j2,]-dat[i,]) )/( sum((dat[j1,]-dat[i,])*(dat[j1,]-dat[i,]))*sum((dat[j2,]-dat[i,])*(dat[j2,]-dat[i,])) )  )
      }
    }
    AS[i]<-var(angles.vect)
  }
  
  AS
}


fastabod.outlier<-function(dat,k=10,standardise=F,verbose=F){
  ## method described in Kriegel et al (2008)
  ## k = number of nearest neighbours

  #library(RANN) # required for nn2 subroutine

  nc<-ncol(dat)
  nr<-nrow(dat)

  if(standardise==TRUE){
    for(i in 1:nc){
      if(sum(is.na(dat[,i]))<nr){dat[,i]<-(dat[,i]-mean(dat[!is.na(dat[,i]),i]))/sd(dat[!is.na(dat[,i]),i])}
    }
  }

  # get list with kNN indices
    nn.list<-nn2(dat=dat,query=dat,k=(k+1),eps=0)

  AS<-rep(NA,nr)
  for(i in 1:nr){
    angles.vect<-numeric(0)
    for(j1 in nn.list$nn.idx[i,-1]){
      for(j2 in (nn.list$nn.idx[i,-1])[nn.list$nn.idx[i,-1]!=j1]){
        if(verbose){print(paste("SAMPLE = ",i,", j1 = ",j1,", j2 = ",j2,sep=""))}
        angles.vect<-c(  angles.vect,sum( (dat[j1,]-dat[i,])*(dat[j2,]-dat[i,]) )/( sum((dat[j1,]-dat[i,])*(dat[j1,]-dat[i,]))*sum((dat[j2,]-dat[i,])*(dat[j2,]-dat[i,])) )  )
      }
    }
    AS[i]<-var(angles.vect)
  }
  AS
}

aggarwal.outlier<-function(dat,phi,k=(-99),thr=0.01,standardise=F,verbose=F){
  ##the method from Aggarwal and Yu (2005)

  nc<-ncol(dat)
  nr<-nrow(dat)

  if(standardise==TRUE){
    for(i in 1:nc){
      if(sum(is.na(dat[,i]))<nr){dat[,i]<-(dat[,i]-mean(dat[!is.na(dat[,i]),i]))/sd(dat[!is.na(dat[,i]),i])}
    }
  }

  if(k==-99){k<-floor(log(nr/10)/log(phi))}
  if(verbose){print(k)}

  miss<-vector("list",nc)
  grid<-vector("list",nc)

  AS<-rep(0,nr)

  ## setting up the grid
  for(j in 1:nc){
    min.tmp<-min(dat[!is.na(dat[,j]),j])
    max.tmp<-max(dat[!is.na(dat[,j]),j])
    nbin<-round(sum(!is.na(dat[,j]))/phi)

    grid[[j]]<-vector("list",2)
    names(grid[[j]])<-c("boundaries","bins")
    grid[[j]]$bins<-vector("list",phi)
    grid[[j]]$boundaries<-min.tmp
    sorted.tmp<-sort(dat[!is.na(dat[,j]),j])
    for(i in 1:(phi-1)){
      grid[[j]]$boundaries<-c(grid[[j]]$boundaries,sorted.tmp[i*nbin])
      grid[[j]]$bins[[i]]<-(1:nr)[!is.na(dat[,j]) & dat[,j]>=grid[[j]]$boundaries[i] & dat[,j]<grid[[j]]$boundaries[i+1]]      
    }
    grid[[j]]$boundaries<-c(grid[[j]]$boundaries,max.tmp)
    grid[[j]]$bins[[phi]]<-(1:nr)[!is.na(dat[,j]) & dat[,j]>=grid[[j]]$boundaries[phi]]      
  }
  AS.tmp2<-numeric(0)
  ## computing the sparsity score for each cube and, if sparse, flag the object in the cube
  subsets<-vector("list")
  check1<-1
  #looping through the C(n,k) subspaces
  while(check1<=as.integer(round(combinatorial(nc,k)))){
    subsets[[check1]]<-sort(sample(x=1:nc,size=k,replace=F))
    if(length(unique(subsets))==check1){
      #looping through the phi^k k-dimensional cubes
	check2<-1
      cubes<-vector("list")
	while(check2<=phi^k){
        cubes[[check2]]<-sample(x=1:phi,size=k,replace=T) 
	  if(length(unique(cubes))==check2){
	    obj_in_cube<-integer(0)
	    for(i in 1:k){
	      obj_in_cube<-c(obj_in_cube,grid[[subsets[[check1]][i]]]$bins[[cubes[[check2]][i]]])
	    }
          obj_in_cube2<-obj_in_cube
          for(i in 1:length(obj_in_cube)){
            if(sum(obj_in_cube==obj_in_cube[i])!=k){obj_in_cube2<-obj_in_cube2[obj_in_cube2!=obj_in_cube[i]]}
	    }
	    obj_in_cube<-unique(obj_in_cube2)
          nobj<-length(obj_in_cube)
          rm(obj_in_cube2)

	    AS.tmp<-(nobj-nr*(1/phi)^k)/sqrt(nr*((1/phi)^k)*(1-(1/phi)^k))
          AS.tmp2<-c(AS.tmp2,AS.tmp)
          if(nobj>0 & AS.tmp<qnorm(thr)){AS[obj_in_cube]<-1}
	    check2<-check2+1;
	  }
	}
      check1<-check1+1; 
    }
  }
  AS

}


combinatorial<-function(n,k){
  if(trunc(n)!=n | trunc(k)!=k | n<0 | k<0 | k>n){stop("n,k need to be non-negative integers, and n>=k")}
  if(n<150){
  	res<-factorial(n)/(factorial(k)*factorial(n-k))
  }
  if(n>=150){
  	prod1<-1
  	for(i in 1:k){prod1<-prod1*(n-k+i)}
  	res<-prod1/factorial(k)
  }
  res	
}


quantile.opt<-function(x,probs,type=7){
	x<-x[!is.na(x) & abs(x)<Inf]
	quantile(x=x,probs=probs,type=type)
}

sample.opt<-function(x,size,replace=F){
  if(replace==F & size>length(x)){res<-sample(x=x,size=size,replace=T)}
  if(replace==F & size<=length(x)){res<-sample(x=x,size=size,replace=F)}
  if(replace==T){res<-sample(x=x,size=size,replace=T)}
  res
}


quantile.emp<-function(x,X,append=1,type=1){
	# computes empirical quantiles
	# x = object of interest (whose empirical quantile is computed) NB x can be a vector itself
	# X = reference vector used to compute the emirical quantile
	# append = 0 (x already contained in X) or 1 (append x to X)
	# type = 1: sum_i( I(X[i]<=x)/n )
	# type = 2: sum_i( I(X[i]<=x)/(n+1) ) NB will always be <1
	# type = 3: sum_i( (I(X[i]<=x)+1)/(n+2) ) NB will always be >0 and <1
	
	if(length(x)==1){
	  if(append==0){X<-X[!is.na(X) & abs(X)<Inf]}
      if(append==1){X<-c(x,X[!is.na(X) & abs(X)<Inf])}
	    if(is.na(x)){res<-NA}
	    if(!is.na(x)){
		    if(type==1){res<-sum(X<=x)/length(X)}
	  	    if(type==2){res<-sum(X<=x)/(length(X)+1)}
		    if(type==3){res<-(sum(X<=x)+1)/(length(X)+2)}
	    }
	}
	if(length(x)>1){
	  res<-rep(NA,length(x))
	  if(append==0){X<-X[!is.na(X) & abs(X)<Inf]}
      if(append==1){X<-c(x,X[!is.na(X) & abs(X)<Inf])}
	  for(i in 1:length(x)){
	      if(!is.na(x[i])){
		      if(type==1){res[i]<-sum(X<=x[i])/length(X)}
	  	      if(type==2){res[i]<-sum(X<=x[i])/(length(X)+1)}
		      if(type==3){res[i]<-(sum(X<=x[i])+1)/(length(X)+2)}
	      }
	  }
	}	
	res
}


smoothBootStrapMedian<-function(smpl,N=100,probs=c(0.5,0.99),noiseSD=(1/sqrt(length(smpl)))){
	res<-matrix(numeric(length=3*N),ncol=3)
	for(i in 1:N){
		BS<-sample(x=smpl,size=length(smpl),replace=T)+rnorm(n=length(smpl),mean=0,sd=noiseSD)
		res[i,1]<-quantile.opt(x=BS,probs=probs[1])
		res[i,2]<-quantile.opt(x=BS,probs=probs[2])
    }
    rm(BS)
    res
}


smoothBootStrapMedianAlt<-function(data,Nboot=100,probs=c(0.5,0.99),noiseSD=-99,N=5,pTail=0.99,noAS){
    # data = matrix with ASs of sources with at least noAS ASs
    # noAS = number of finite, non-missing ASs each object is required to have
    # noiseSD = std of Gaussian noise in smooth bootstrap sample; if <0 it will be set to 1/sqrt(n_samp)
    # Nboot = number of bootstrap samples
    # N = N as in TopN AS combination scheme
    # pTail = pTail as in TopTail AS combination scheme
    # probs = quantiles of interest (0.5 = median etc)
      res<-vector("list",3)
      names(res)<-c("AD.Ext","AD.Ave","AD.TopTail")
	for(i in 1:3){res[[i]]<-matrix(numeric(length=3*Nboot),ncol=3)}
	for(i in 1:Nboot){
            smpl<-dataprep.bootstrap(data,noAS=noAS,pTail=pTail,N=N)
            if(noiseSD<0){noiseSTD<-1/(sqrt(length(smpl)))}
            if(noiseSD>=0){noiseSTD<-noiseSD}
		BS<-sample(x=smpl$AD.Ext,size=length(smpl),replace=T)+rnorm(n=length(smpl),mean=0,sd=noiseSTD)
		res$AD.Ext[i,1]<-quantile.opt(x=BS,probs=probs[1])
		res$AD.Ext[i,2]<-quantile.opt(x=BS,probs=probs[2])
		
		sample.median<-res$AD.Ext[i,1]
            dens.emp<-density(BS[!is.na(BS) & abs(BS)<Inf],kernel="gaussian",from=sample.median,to=sample.median,n=1)$y
              # this is the empirical density function evaluted at the sample median
            # NB since p=1-p=0.5, all the influence functions for the median statistic are the same, namely 0.5/dens.emp(sample_median)
            infl.fun<-0.5/dens.emp
            res$AD.Ext[i,3]<-(infl.fun^2)/nrow(CMsamp)
		
            rm(BS,sample.median,infl.fun,dens.emp)
		BS<-sample(x=smpl$AD.Ave,size=length(smpl),replace=T)+rnorm(n=length(smpl),mean=0,sd=noiseSTD)
		res$AD.Ave[i,1]<-quantile.opt(x=BS,probs=probs[1])
		res$AD.Ave[i,2]<-quantile.opt(x=BS,probs=probs[2])
		
		sample.median<-res$AD.Ave[i,1]
            dens.emp<-density(BS[!is.na(BS) & abs(BS)<Inf],kernel="gaussian",from=sample.median,to=sample.median,n=1)$y
              # this is the empirical density function evaluted at the sample median
            # NB since p=1-p=0.5, all the influence functions for the median statistic are the same, namely 0.5/dens.emp(sample_median)
            infl.fun<-0.5/dens.emp
            res$AD.Ave[i,3]<-(infl.fun^2)/nrow(CMsamp)
		
            rm(BS,sample.median,infl.fun,dens.emp)
		BS<-sample(x=smpl$AD.TopTail,size=length(smpl),replace=T)+rnorm(n=length(smpl),mean=0,sd=noiseSTD)
		res$AD.TopTail[i,1]<-quantile.opt(x=BS,probs=probs[1])
		res$AD.TopTail[i,2]<-quantile.opt(x=BS,probs=probs[2])
		
		sample.median<-res$AD.TopTail[i,1]
            dens.emp<-density(BS[!is.na(BS) & abs(BS)<Inf],kernel="gaussian",from=sample.median,to=sample.median,n=1)$y
              # this is the empirical density function evaluted at the sample median
            # NB since p=1-p=0.5, all the influence functions for the median statistic are the same, namely 0.5/dens.emp(sample_median)
            infl.fun<-0.5/dens.emp
            res$AD.TopTail[i,3]<-(infl.fun^2)/nrow(CMsamp)
		
            rm(BS,sample.median,infl.fun,dens.emp)
    }
    res
}


dataprep.bootstrap<-function(data,noAS,combfun,...){
	  # combines a matrix with AS to a single AS vector using combination function combfun, but with only j As per object
        # data = matrix with ASs of sources with at least noAS ASs
        # noAS = number of finite, non-missing ASs each object is required to have
    noAS<-as.integer(round(noAS))
    data.smpl<-rep(NA,length=noAS*nrow(data))
    data.smpl<-matrix(data.smpl,nrow=nrow(data))
    
    for(i in 1:nrow(data)){
      x.smpl<-(1:ncol(data))[!is.na(data[i,]) & abs(data[i,])<Inf]
      if(sum(!is.na(data[i,]))<noAS){data.smpl[i,]<-rep(NA,length=noAS)}
      if(sum(!is.na(data[i,]))>=noAS){data.smpl[i,]<-data[i,sample.opt(x=x.smpl,size=noAS,replace=F)]}
      rm(x.smpl)
    }
    rm(data)
    res<-AnomDet.combine(data.smpl,combfun=combfun,...)
    res
}


varMV<-function(dat){
  nr<-nrow(dat)
  nc<-ncol(dat)

  res<-matrix(rep(NA,length=(nc*nc)),ncol=nc)

  for(i in 1:nc){
    if(i>1){
      for(j in 1:(i-1)){
        res[i,j]<-cov(dat[!is.na(dat[,i]) & !is.na(dat[,j]),i],dat[!is.na(dat[,i]) & !is.na(dat[,j]),j])
        res[j,i]<-res[i,j]
      }
    }
    res[i,i]<-var(dat[!is.na(dat[,i]),i])
  }

  res

}


hist.char<-function(x,cats=NULL){
  ## computes counts for character vectors which can then be used to plot a histogram
  ## x = character vector
  ## cats = optional argument listing the categories to check for 

  if(!is.character(x)){stop("Input needs to be of type character.")}

  if(is.null(cats)){cats<-sort(unique(x))}
  counts<-integer(length=length(cats))

  for(i in 1:length(cats)){
    counts[i]<-sum(x==cats[i])
  }

  res<-list(cats,counts)
  names(res)<-c("categories","counts")

  res
}











