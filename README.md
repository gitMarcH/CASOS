# CASOS

This is the R implementation of the CASOS algorithm (https://doi.org/10.1002/sam.11167).

The algorithm was (and still is) previously available from my personal website and is old (2012) R code which I have decided to migrate to GitHub for archiving.
The code should still run with recent (4.1.0) versions of R, but no guarantee -- I have not tried it.

Original README copied below.

## INSTALLATION
To load CASOS, simply source the file "CASOSsource.R" from within R:
source("/path/CASOSsource.R")
where "path" is the directory in which you have saved the R file.

You will need to install the RANN library (needed for finding the nearest neighbours for LOF).


## USAGE:

To run CASOS is a two-stage process: first you need to run CASOS(...) to compute anomaly scores in all subspaces of a certain dimensionality, then CASOS.combine(...) to combine these to a single anomaly score per object.
If LOF is used as anomaly score computation method, then the AS matrix for subspaces of dimensionality np (parameter D from the CASOS paper) is the np^th element of the list $LOF.sc, and can be retrieved by $LOF.sc[[np]].

Example (for np (i.e. D) = 3):

ASresults<-CASOS(X=data,np=3,k=100)
AScombined<-CASOS.combine(ASresults$LOF.sc[[3]],combfun=Comb.Ext)


## PARAMETERS:

  (a) CASOS(X,X.weight=NULL,weight.scheme=NULL,np=3,ASmethod="LOF",k=5,h=1,standardise=TRUE,do.all.dim=FALSE)

  You will mostly just need to set your input data (X), the number of nearest neighbours (k) and the maximum number of dimensions of subspaces (np).
  There are more parameters that can be set, they are listed below.

	# X = data matrix;
	#     all columns should be numerical variables
	# X.weight = matrix with weights per variable
	#            will be used to compute the weights to apply to each subspace
	#            if NULL then all weights are equal and set to 1
	#	(defaults to NULL)
	# weight.scheme = how to use X.weight to weight the anomaly score (one of "ratio1","ratio2","ratio0-5","ratioAVG", "exp1","exp2","exp0-5" cf. comments in code for more info)
	#	(defaults to NULL)
	# np = maximum dimensionality of projections;
	#      all projections of all dimensinalities<=np will be computed if do.al.dim==T, otherwise only projections of dimensionality exactly np will be computed
	# ASmethod = method to be used to compute anomaly scores;
	#            for the moment are supported:
	#               "LOF" (Local Outlier Factor; Breuning et al. ACM SIGMOD Conf. Proc.; 2000)
	#               "LDF" (Outlier Detection with Kernel Density Functions; Latecki, Lazarevic and Pokrajac, LNAI 4571; 2007)
	#               "LDF.ecld" (a version of LDF using only Euclidean distance rather than Euclidean to find the kNN and then Mahalanobis for score calculation)
	#	(defaults to "LOF")
	# k = number of nearest neighbours (for ASmethod="LOF", "LDF" and "LDF.ecld")
	# h = smoothing parameter (for ASmethod="LDF" and "LDF.ecld")
	# standardise = logical;
	#               if TRUE all columns are standardised to zero mean and unit variance
	#	(defaults to TRUE)
	# do.all.dim = logical;
	#         if TRUE then AS in all subspaces with dim<=np will be computed, else only subspaces of dimesnionality np will be computed
	#	(defaults to FALSE)

  (B) CASOS.combine(ASmat,combfun,...)

  You will need to set the matrix of anomaly score (ASmat) and the combination function to be used (combun). If the combination function requires further input parameters, these need to be set as well.

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
