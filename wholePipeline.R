
mvslPipeline<-function(disease,treefile,varnames,Merrorvarnames,vtorun=rep(TRUE,4),numSims=200000,directCompare=FALSE,IndepDiff=FALSE,simset=1,pvalcutoff=0.1,maxnumsigs=2){
    controldatafile<-paste(disease,"Normal.csv",sep="")
    casedatafile<-paste(disease,"Patient.csv",sep="")

    if (IndepDiff){directCompare<-FALSE}
    outprefix<-paste(disease,"Normal",sep="")

    modelfile<-paste(outprefix,"_BestModel.RData",sep="")
    if (vtorun[1]){FindBestModel(treefile,controldatafile,varnames,Merrorvarnames,outprefix)}
	load(modelfile)
        modelParams<-BestModel$BestModel$ParamsInModel
	evolmodel<-"BM"
        if (!is.null(modelParams$A)){evolmodel<-"OU"}
	if (!is.null(modelParams$B)){evolmodel<-"SLOUCH"}
    simulateddatafile<-paste(outprefix,"_",evolmodel,"sim_",simset,".RData",sep="")

    if (vtorun[2]){SimulateHealthyModel(modelfile,treefile,simulateddatafile,numSims)}
    nulldistributionfile<-NA
    if (directCompare){
	nulldistributionfile<-paste(outprefix,"_",evolmodel,"NullDistributionDirectCompare_",simset,".RData",sep="")
    }else{
        if (IndepDiff){nulldistributionfile<-paste(outprefix,"_",evolmodel,"NullDistributionIndepDiff_",simset,".RData",sep="")}
	else{nulldistributionfile<-paste(outprefix,"_",evolmodel,"NullDistributionConditionalHealthy_",simset,".RData",sep="")}
    }

    if (vtorun[3]){getNullDistribution(controldatafile,casedatafile,treefile,simulateddatafile,varnames,Merrorvarnames,nulldistributionfile,directCompare,IndepDiff)}
    significanttipsfile<-NA
    if (directCompare){
	significanttipsfile<-paste(outprefix,"_",evolmodel,"SignificantReceptorsDirectCompare_",simset,".RData",sep="")
    }else{
	if (IndepDiff){significanttipsfile<-paste(outprefix,"_",evolmodel,"SignificantReceptorsIndepDiff_",simset,".RData",sep="")}
	else{significanttipsfile<-paste(outprefix,"_",evolmodel,"SignificantReceptorsConditionalHealthy_",simset,".RData",sep="")}
    }

    if (vtorun[4]){CompareToNull(controldatafile,casedatafile,nulldistributionfile,varnames,significanttipsfile,directCompare)}

    signiffile<-paste(disease,"SignificantReceptors",sep="")
    if (directCompare){
	signiffile<-paste(signiffile,"DirectCompare",sep="")    
    }else{
	if (IndepDiff){signiffile<-paste(signiffile,"IndepDiff",sep="")}
	else{signiffile<-paste(signiffile,"ConditionalHealthy",sep="")}
    }
    signiffile<-paste(signiffile,".txt",sep="")    
    if (vtorun[5]){
	sink(signiffile);signifres<-FindSignif(disease,significanttipsfile,pvalcutoff);sink()
	if (vtorun[6]){## Now find groups
	    #signiffilegroups<-paste(signiffile,"Groups.txt",sep="")    
	    sink(signiffile,append=TRUE)
	    cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	    cat("\n");	    cat("\n");
	    cat("Significantly different groups")
	    cat("\n")
	    significanttipgroupsfile<-NA
	    if (directCompare){
		significanttipgroupsfile<-paste(outprefix,"_",evolmodel,"SignificantReceptorsGroupDirectCompare_",simset,".RData",sep="")
	    }else{
		if (IndepDiff){significanttipgroupsfile<-paste(outprefix,"_",evolmodel,"SignificantReceptorsGroupIndepDiff_",simset,".RData",sep="")}
		else{significanttipgroupsfile<-paste(outprefix,"_",evolmodel,"SignificantReceptorsGroupConditionalHealthy_",simset,".RData",sep="")}
	    }
	    CompareToNullGroups(controldatafile,casedatafile,nulldistributionfile,varnames,signifres$SignificantFound,significanttipgroupsfile,directCompare,maxnumsigs=maxnumsigs)
	    PrintSignifGroups(disease,significanttipgroupsfile,pvalcutoff)
	    sink()
	}
    }
}


