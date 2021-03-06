evol.model.est<-function(phyltreeouch,dfdata,regimes=NULL,root.regime=NULL,M.error=NULL,repeats=3,model.setups=NULL,predictors=NULL,kY=NULL,doPrint=FALSE){
    BestModel<-list(BestModel=NA,aic.c=10000,i=NA,model=NA,evolmodel=NA)
    testedModels<-list()
    j<-1
    if (is.null(model.setups)){
## List below defines all possible interesting models
	model.setups<-generate.model.setups()    
    }
    for (i in 1:repeats){
	for (k in 1:length(model.setups)){
	    if (model.setups[[k]]$evolmodel=="BM"){
		if (doPrint){print("Doing estimation for BM model")}
		BMres<-NULL
		tryCatch({
		    BMres<-BrownianMotionModel(phyltree=phyltreeouch,data=dfdata,predictors=predictors,M.error=M.error,calcCI=FALSE)
		},error=function(e){print(e)})
		testedModels[[j]]<-list()
		testedModels[[j]]$result<-BMres;testedModels[[j]]$aic.c<-NA;testedModels[[j]]$model<-model.setups[[k]];j<-j+1	    
		if (!(is.null(BMres))&&(BMres$ParamSummary$aic.c < BestModel$aic.c)){
		    BestModel$BestModel<-BMres
		    BestModel$aic.c<-BMres$ParamSummary$aic.c
		    BestModel$i<-i
		    BestModel$model.call<-"BM"
		    BestModel$model<-model.setups[[k]]
		    BestModel$evolmodel<-"BM"
		}
	    }
	    if (model.setups[[k]]$evolmodel=="ouch"){
		if(doPrint){print(paste("Doing estimation for ouch model with A: ",model.setups[[k]]$Atype," with diagonal: ",model.setups[[k]]$diagA," Syy: ",model.setups[[k]]$Syytype,sep=""))}
		OUres<-NULL
		tryCatch({
		    OUres<-ouchModel(phyltree=phyltreeouch,data=dfdata,regimes=regimes,regimes.times=NULL,root.regime=root.regime,predictors=predictors,M.error=M.error,Atype=model.setups[[k]]$Atype,Syytype=model.setups[[k]]$Syytype,calcCI=FALSE,diagA=model.setups[[k]]$diagA)
		},error=function(e){print(e)})
		testedModels[[j]]<-list()
		testedModels[[j]]$result<-OUres;testedModels[[j]]$aic.c<-NA;testedModels[[j]]$model<-model.setups[[k]];j<-j+1
    		if (!is.null(OUres)){
    		    if (is.list(OUres$MaxLikFound)){
			if (OUres$MaxLikFound$ParamSummary$aic.c < BestModel$aic.c){
			    BestModel$BestModel<-OUres$MaxLikFound
			    BestModel$aic.c<-OUres$MaxLikFound$ParamSummary$aic.c
			    BestModel$i<-i
			    BestModel$evolmodel<-"ouch"
			    BestModel$model.call<-paste("ouch model with A: ",model.setups[[k]]$Atype," with diagonal: ",model.setups[[k]]$diagA," Syy: ",model.setups[[k]]$Syytype,sep="")
			    BestModel$model<-model.setups[[k]]
			}
		    }else{
			if (OUres$FinalFound$ParamSummar$aic.c < BestModel$aic.c){
			    BestModel$BestModel<-OUres$FinalFound
			    BestModel$aic.c<-OUres$FinalFound$ParamSummar$aic.c
			    BestModel$i<-i
			    BestModel$evolmodel<-"ouch"
			    BestModel$model.call<-paste("ouch model with A: ",model.setups[[k]]$Atype," with diagonal: ",model.setups[[k]]$diagA," Syy: ",model.setups[[k]]$Syytype,sep="")
			    BestModel$model<-model.setups[[k]]
			}
		    }
		}
	    }
	    if (model.setups[[k]]$evolmodel=="mvslouch"){
		if (is.null(kY)){
		    dfdata.mvsl<-dfdata[,c(setdiff(1:ncol(dfdata),predictors),predictors)]
		    kY<-ncol(dfdata)-predictors
		}else{dfdata.mvsl<-dfdata}
		if(doPrint){print(paste("Doing estimation for mvlslouch model with A: ",model.setups[[k]]$Atype," with diagonal: ",model.setups[[k]]$diagA," Syy: ",model.setups[[k]]$Syytype,sep=""))}
	    	mvslres<-NULL
	    	tryCatch({
	    	    mvslres<-mvslouchModel(phyltree=phyltreeouch,data=dfdata.mvsl,kY,regimes=regimes,regimes.times=NULL,root.regime=root.regime,predictors=predictors,M.error=M.error,Atype=model.setups[[k]]$Atype,Syytype=model.setups[[k]]$Syytype,calcCI=FALSE,diagA=model.setups[[k]]$diagA)
		},error=function(e){print(e)})
		testedModels[[j]]<-list()
		testedModels[[j]]$result<-mvslres;testedModels[[j]]$aic.c<-NA;testedModels[[j]]$model<-model.setups[[k]];j<-j+1
    		if (!is.null(mvslres)){
		    if (is.list(mvslres$MaxLikFound)){
			if (mvslres$MaxLikFound$ParamSummar$aic.c < BestModel$aic.c){
			    BestModel$BestModel<-mvslres$MaxLikFound
			    BestModel$aic.c<-mvslres$MaxLikFound$ParamSummary$aic.c
			    BestModel$i<-i
			    BestModel$model.call<-paste("mvslouch model with A: ",model.setups[[k]]$Atype," with diagonal: ",model.setups[[k]]$diagA," Syy: ",model.setups[[k]]$Syytype,sep="")
			    BestModel$model<-model.setups[[k]]
			    BestModel$evolmodel<-"mvslouch"
			}
		    }else{
			if (mvslres$FinalFound$ParamSummar$aic.c < BestModel$aic.c){
			    BestModel$BestModel<-mvslres$FinalFound
			    BestModel$aic.c<-mvslres$FinalFound$ParamSummary$aic.c
			    BestModel$i<-i
			    BestModel$model.call<-paste("mvslouch model with A: ",model.setups[[k]]$Atype," with diagonal: ",model.setups[[k]]$diagA," Syy: ",model.setups[[k]]$Syytype,sep="")
			    BestModel$model<-model.setups[[k]]
			    BestModel$evolmodel<-"mvslouch"
			}
		    }
		}
	    }
	}
    }   
    ballPosEig<-FALSE 
    numtraits<-1
    if ((BestModel$evolmodel=="ouch")||(BestModel$evolmodel=="mvslouch")){
	if (length(which(Re(BestModel$BestModel$ParamSummary$phyl.halflife$halflives["halflife",])>0)) == length(BestModel$BestModel$ParamSummary$phyl.halflife$halflives["halflife",])){ballPosEig<-TRUE}
	numtraits<-ncol(BestModel$BestModel$ParamsInModel$A)
    }
    if (BestModel$evolmodel=="BM"){numtraits<-ncol(BestModel$BestModel$ParamsInModel$Sxx)}    
    BestModel$model.description<-.generate.model.description(BestModel$evolmodel,BestModel$model,numtraits,ballPosEig)
    BestModel$key.properties<-.extract.model.key.properties(BestModel$evolmodel,BestModel$BestModel,k=numtraits)
    BestModel$parameter.SE<-.generate.model.se(BestModel$evolmodel)
    BestModel<-BestModel[c("model.description","key.properties","parameter.SE","aic.c","evolmodel","model","model.call","BestModel","i")]
    BestModel$BestModel<-BestModel$BestModel[c("ParamsInModel","ParamSummary","LogLik","HeuristicSearchPointFinalFind")]
    list(BestModel=BestModel,testedModels=testedModels,model.setups=model.setups,repeats=repeats)
}


generate.model.setups_2<-function(){
	list(
	    list(evolmodel="BM"),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="Diagonal",diagA="Positive")
	)
}
	

generate.model.setups<-function(){
## This function should be really hidden but is made available so that the user can create a personalized version of it to speed up estimation.
	list(
	    list(evolmodel="BM"),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="SymmetricPositiveDefinite",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="UpperTri",diagA="NULL"),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="SymmetricPositiveDefinite",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="Diagonal",diagA="NULL"),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="SymmetricPositiveDefinite",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA="NULL"),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SymmetricPositiveDefinite",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="UpperTri",diagA="NULL"),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SymmetricPositiveDefinite",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="Diagonal",diagA="NULL"),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SymmetricPositiveDefinite",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA="NULL"),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA="Positive")
	)
}

.generate.model.description<-function(evolmodel,model.setup,numtraits,ballPosEig){
    chdesc<-"This is a qualitative description of the best found model by the estimation package. See also the field $key.properties for quantative descriptions. "
    if (evolmodel=="BM"){	
	chdesc<-paste(chdesc,"Brownian motion model, no stabilizing selection observed. ",sep="")
	if (numtraits==1){chdesc<-paste(chdesc,"Phenotype randomly oscillates ",sep="")}
	else{chdesc<-paste(chdesc,"Phenotypes randomly oscillate ",sep="")}
	chdesc<-paste(chdesc,"around the initial root state accumulating more and more variance with time.",sep="")
    }
    if ((evolmodel=="ouch") || (evolmodel=="mvslouch")){    

	if (numtraits==1){
	    if ((evolmodel=="ouch") && ballPosEig){chdesc<-paste(chdesc,"Single trait under stabilizing selection, adapting to a deterministic optimum. ",sep="")}
	    if ((evolmodel=="ouch") && !ballPosEig){chdesc<-paste(chdesc,"Single trait diverging. Results needs careful interpretation ",sep="")}
	    if ((evolmodel=="mvslouch") && ballPosEig){chdesc<-paste(chdesc,"Single trait under stabilizing selection, adapting to a randomly evolving optimum. ",sep="")}
	    if ((evolmodel=="mvslouch") && !ballPosEig){chdesc<-paste(chdesc,"Single trait diverging from a randomly evolving value. Results needs careful interpretation ",sep="")}
	}else{
	    if (!ballPosEig){
		if (evolmodel=="ouch") {chdesc<-paste(chdesc,"Multiple traits diverging. Results needs careful interpretation ",sep="")}
		if (evolmodel=="mvslouch") {chdesc<-paste(chdesc,"Multiple traits diverging from a randomly evolving value. Results needs careful interpretation ",sep="")}
	    }else{
		if (evolmodel=="ouch") {chdesc<-paste(chdesc,"Multiple traits under stabilizing selection, adapting to a deterministic optimum. ",sep="")}
		if (evolmodel=="mvslouch") {chdesc<-paste(chdesc,"Multiple traits under stabilizing selection, adapting to s randomly evolving optimum. ",sep="")}
		if ((model.setup$Atype== "Diagonal")|| (model.setup$Atype=="SingleValueDiagonal")){
		    if ((model.setup$Syytype== "Diagonal")|| (model.setup$Syytype=="SingleValueDiagonal")){		
			chdesc<-paste(chdesc,"All traits are evolving independently. They do not effect each others' adaptation to the primary optimum nor are their random perturbations correlated. ",sep="")
			if (model.setup$Atype=="SingleValueDiagonal"){chdesc<-paste(chdesc,"All traits have the same rate of adaptation to the primary optimum. ",sep="")}
			if (model.setup$Syytype=="SingleValueDiagonal"){chdesc<-paste(chdesc,"All traits experiance the same magnitude of random perturbations. ",sep="")}		    
		    }else{
			chdesc<-paste(chdesc,"All traits are adapting independently to the primary optimum. All dependencies between traits are due to the random perturbations. ",sep="")
			if (model.setup$Atype=="SingleValueDiagonal"){chdesc<-paste(chdesc,"All traits have the same rate of adaptation to the primary optimum. ",sep="")}		    
		    }		
		}
		else{
		    if ((model.setup$Atype== "UpperTri")|| (model.setup$Atype=="LowerTri")){
			chdesc<-paste(chdesc,"We have a clear causality pattern in traits' adaptation and primary optimum. ",sep="")
			if (model.setup$Atype== "UpperTri"){
			    chdesc<-paste(chdesc,"The ''bottom'' trait in the A matrix is effecting the primary optimum of all other traits. The second last one all the traits' except the bottom one's and so on. The ''top'' trait's primary optimum is effected by all other traits. ",sep="")
			}else{
			    chdesc<-paste(chdesc,"The ''top'' trait in the A matrix is effecting the primary optimum of all other traits. The second one all the traits' except the top one's and so on. The ''bottom'' trait's primary optimum is effected by all other traits. ",sep="")
			}
		    }else{
			chdesc<-paste(chdesc,"There is no clear causality pattern in traits' adaptation. All traits seem to effect each other's primary optimum. However look in the A matrix for 0 or very small values. These will indicate bivariate causality effects. ",sep="")
		    }
		    if ((model.setup$Syytype== "Diagonal")|| (model.setup$Syytype=="SingleValueDiagonal")){		
			chdesc<-paste(chdesc,"The random perturbations are not correlated between different traits. All dependencies between the traits come from their co-adaptation. ",sep="")
			if (model.setup$Syytype=="SingleValueDiagonal"){chdesc<-paste(chdesc,"All traits experiance the same magnitude of random perturbations. ",sep="")}		    
		    }else{
			chdesc<-paste(chdesc,"The random perturbations are correlated between different traits. The traits are dependent both through their co-adaptation and correlated random perturbations. ",sep="")		    
		    }
		}
	    }
	}
    }    
    chdesc
}

.generate.model.se<-function(evolmodel){
    chsedesc<-"The function did not calculate any confidence intervals as this takes a long time. You have to do it manually. Here is code to call: "
    if (evolmodel=="BM"){
        chsedesc<-paste(chsedesc,"SummarizeBM(phyltree=your.phylogenetic.tree,data=your.data.frame,modelParams=your.estimated.model.parameters,dof=your.degrees.of.freedom,M.error=your.measurement.error,predictors=your.predictors,calcCI=TRUE) ",sep="")
    }
    if (evolmodel=="ouch"){
	chsedesc<-paste(chsedesc,"SummarizeOUCH(phyltree=your.phylogenetic.tree,data=your.data.frame,modelParams=your.estimated.model.parameters,regimes=your.regimes,dof=your.degrees.of.freedom,M.error=your.measurement.error,predictors=your.predictors,Atype=your.Atype,Syytype=your.Syytype,calcCI=TRUE) ", sep="")
    }
    if (evolmodel=="mvslouch"){
	chsedesc<-paste(chsedesc,"SummarizeMVSLOUCH(phyltree=your.phylogenetic.tree,data=your.data.frame,modelParams=your.estimated.model.parameters,regimes=your.regimes,dof=your.degrees.of.freedom,M.error=your.measurement.error,predictors=your.predictors,Atype=your.Atype,Syytype=your.Syytype,calcCI=TRUE) ", sep="")
    }
    chsedesc<-paste(chsedesc,"where your.phylogenetic.tree is your phylogeny in ouch format (passed to the parameter phyltreeouch when calling evol.model.est), your.data.frame is your trait data (passed to the parameter dfdata when calling evol.model.est), your.estimated.model.parameters are the estimated parameters (can be found in the slot $BestModel$BestModel$ParamsInModel in the returned object),  your.degrees.of.freedom are the degrees of freedom of the best found model (can be found in the slot $BestModel$BestModel$ParamSummary$dof in the returned object), your.measurement.error is your measurement error structure (passed to the parameter M.error when calling evol.model.est, may be left NULL), your.predictors are your indicated predictor variables (passed to the parameter predictors when calling evol.model.est, may be left NULL)",sep="")    
    if ((evolmodel=="ouch") || (evolmodel=="mvslouch")){
	chsedesc<-paste(chsedesc,", your.regimes is your regimes vector (passed to the parameter regimes when calling evol.model.est), your.Atype is the class to which the A matrix belongs in the best found model (can be found in the slot $BestModel$model$Atype in the returned object), your.Syytype is the class to which the A matrix belongs in the best found model (can be found in the slot $BestModel$model$Syytype in the returned object). ", sep="")
    }else{chsedesc<-paste(chsedesc,".",sep="")}
    chsedesc
}

.extract.model.key.properties<-function(evolmodel,estimated.model,k=NULL){
    lkey.properties<-list()    
    if (evolmodel=="BM"){lkey.properties$evolution.model<-"Brownian motion model"}
    if ((evolmodel=="ouch")||(evolmodel=="mvslouch")){
	if ((evolmodel=="ouch") && (k==1)){lkey.properties$evolution.model<-"Ornstein-Uhlenbeck model"}
	if ((evolmodel=="ouch") && (k>1)){lkey.properties$evolution.model<-"Ornstein-Uhlenbeck-Ornstein-Uhlenbeck model"}
	if (evolmodel=="mvslouch"){lkey.properties$evolution.model<-"Ornstein-Uhlenbeck-Brownian motion model"}
	lkey.properties$phylogenenetic.halflives<-list()
	lkey.properties$phylogenenetic.halflives$halflives<-estimated.model$ParamSummary$phyl.halflife$halflives["halflife",]
	lkey.properties$phylogenenetic.halflives$comment<-"Phylogenetic half-lives in the eigenvector (NOT trait) directions. "
	if (length(which(Re(estimated.model$ParamSummary$phyl.halflife$halflives["halflife",])>0))== length(estimated.model$ParamSummary$phyl.halflife$halflives["halflife",])){
	    lkey.properties$phylogenenetic.halflives$comment<-paste(lkey.properties$phylogenenetic.halflives$comment,"All halflives are positive so all traits are under stabilizing selection towards the optimum. See Bartoszek et. al. (2012) for details.",sep="")
	}else{
	    lkey.properties$phylogenenetic.halflives$comment<-paste(lkey.properties$phylogenenetic.halflives$comment,"There are negative halflives so all traits may be diverging. One needs to interpret the results very carefully. See Bartoszek et. al. (2012) for details.",sep="")
	}
	if (!is.null(estimated.model$ParamSummary$evolutionary.regression)){
	    lkey.properties$evolutionary.regression<-list()
	    lkey.properties$evolutionary.regression$regression.coefficents<-estimated.model$ParamSummary$evolutionary.regression
	    lkey.properties$evolutionary.regression$comment<-"The evolutionary regression between traits and the primary optimum for the suite of traits. This is the currently observed relationship. See Bartoszek et. al. (2012) and the help for the mvSLOUCH::mvslouchModel and mvSLOUCH::ouchModel functions for details."
	}
	if (!is.null(estimated.model$ParamSummary$optimal.regression)){
	    lkey.properties$optimal.regression<-list()
	    lkey.properties$optimal.regression$regression.coefficents<-estimated.model$ParamSummary$optimal.regression
	    lkey.properties$optimal.regression$comment<-"The optimal regression between traits and the primary optimum for the suite of traits. This is the currently optimal relationship that the traits want to attain. This however is not always possible, see the discussion in Bartoszek et. al.(2012)."
	}
	if (!is.null(estimated.model$ParamSummary$trait.regression)){
	    lkey.properties$trait.regression<-list()
	    lkey.properties$trait.regression$regression.coefficents<-estimated.model$ParamSummary$trait.regression
	    lkey.properties$trait.regression$comment<-"The regression of each trait on all of the others. See Bartoszek et. al. (2012) and the help for the mvSLOUCH::mvslouchModel and mvSLOUCH::ouchModel functions for details."
	}
	if (!is.null(estimated.model$ParamSummary$corr.matrix)){
	    lkey.properties$corr.matrix<-list()
	    lkey.properties$corr.matrix$correlation.matrix<-estimated.model$ParamSummary$corr.matrix
	    lkey.properties$corr.matrix$comment<-"The currently observed correlation between the traits. See Bartoszek et. al. (2012) and the help for the mvSLOUCH::mvslouchModel and mvSLOUCH::ouchModel functions for details."
	}
	if (!is.null(estimated.model$ParamSummary$stationary.corr.matrix)){
	    lkey.properties$stationary.corr.matrix<-list()
	    lkey.properties$stationary.corr.matrix$correlation.matrix<-estimated.model$ParamSummary$stationary.corr.matrix
	    lkey.properties$stationary.corr.matrix$comment<-"The limiting/stationary correlation between the traits. See Bartoszek et. al. (2012) and the help for the mvSLOUCH::mvslouchModel and mvSLOUCH::ouchModel functions for details."
	}
	lkey.properties$R2<-estimated.model$ParamSummary$RSS$R2
    }
    lkey.properties$aic.c<- estimated.model$BestModel$ParamSummary$aic.c    
    lkey.properties$comment<-"See the slot $BestModel$BestModel$ParamSummary in the returned object for various other hopefully helpful summary statistics and information criteria"
    lkey.properties
}
