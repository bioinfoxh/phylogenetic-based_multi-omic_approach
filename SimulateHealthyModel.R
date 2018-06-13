SimulateHealthyModel<-function(modelfile,treefile,outfile,numSims=200000,M.error=NULL){

    ## preparing the data for the mvSLOUCH package
    phyltreeape<-read.nexus(treefile)
    phyltreeouch<-ape2ouch(phyltreeape,scale=1)
    phyltreeouch@nodelabels[1:(phyltreeouch@nnodes-phyltreeouch@nterm)]<-as.character(1:(phyltreeouch@nnodes-phyltreeouch@nterm))

    load(modelfile)
    modelParams<-BestModel$BestModel$ParamsInModel
    evolmodel<-"bm"
    if (!is.null(modelParams$A)){evolmodel<-"ouch";
	    if(max(diag(modelParams$A))>10){
		diag(modelParams$A)[which(diag(modelParams$A)>10)]<-10    
	    }
	    if(min(diag(modelParams$Syy))<10E-12){
		diag(modelParams$Syy)[which(diag(modelParams$A)<10E-12)]<-1E-12    
	    }
    }
    if (!is.null(modelParams$B)){evolmodel<-"mvslouch"}

    lSimulatedHealthyData<-vector("list",numSims)
    for (i in 1:numSims){
	if (i%%100==1){print(paste("Doing simulation",i))}
	switch(evolmodel,
	    bm={lSimulatedHealthyData[[i]]<-simulBMProcPhylTree(phyltreeouch, modelParams$vX0, modelParams$Sxx, dropInternal = FALSE, M.error=M.error)},
	    ouch={lSimulatedHealthyData[[i]]<-simulOUCHProcPhylTree(phyltreeouch, modelParams, regimes = NULL, regimes.times = NULL, dropInternal = FALSE, M.error=M.error)},
	    mvslouch={lSimulatedHealthyData[[i]]<-simulMVSLOUCHProcPhylTree(phyltreeouch, modelParams, regimes = NULL, regimes.times = NULL, dropInternal = FALSE, M.error=M.error)},
	    {print("The code cannot handle this sort of model")})
    }
    save(lSimulatedHealthyData,phyltreeouch,modelParams,numSims,file=outfile)
}
