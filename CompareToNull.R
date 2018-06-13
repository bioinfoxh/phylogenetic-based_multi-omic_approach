CompareToNull<-function(controldatafile,casedatafile,nulldistributionfile,varnames,outfile,directCompare=FALSE){
    Controldata<-read.csv(controldatafile,header=TRUE,row.names=1)
    Casedata<-read.csv(casedatafile,header=TRUE,row.names=1)
    load(nulldistributionfile) ## here we force the variable name to be lNullRecepDiff
    if (directCompare){StatCaseControl<-Casedata[,varnames]}
    else{StatCaseControl<-(Casedata[,varnames]-Controldata[,varnames])}
    if (length(varnames)==1){
	StatCaseControl<-matrix(StatCaseControl,ncol=1,nrow=length(StatCaseControl))
	colnames(StatCaseControl)<-varnames
    }

    DiffpvalRightTail<-matrix(NA,ncol=ncol(StatCaseControl),nrow=nrow(StatCaseControl))
    DiffpvalLeftTail<-matrix(NA,ncol=ncol(StatCaseControl),nrow=nrow(StatCaseControl))
    DiffpvalAbs<-matrix(NA,ncol=ncol(StatCaseControl),nrow=nrow(StatCaseControl))
    
    for (i in 1:nrow(StatCaseControl)){
	for (j in 1:ncol(StatCaseControl)){
	    if ((!is.na(StatCaseControl[i,j]))&&(!is.na(Casedata[i,j]))&&(!is.na(Controldata[i,j]))){
		irecep<-which(names(lNullRecepDiff)==rownames(StatCaseControl)[i])
		jtrait<-which(names(lNullRecepDiff[[irecep]])==colnames(StatCaseControl)[j])
		DiffpvalAbs[i,j]<-length(which(abs(lNullRecepDiff[[irecep]][[jtrait]])>=StatCaseControl[i,j]))/length(lNullRecepDiff[[irecep]][[jtrait]])
		DiffpvalRightTail[i,j]<-length(which((lNullRecepDiff[[irecep]][[jtrait]])>=StatCaseControl[i,j]))/length(lNullRecepDiff[[irecep]][[jtrait]])
		DiffpvalLeftTail[i,j]<-length(which((lNullRecepDiff[[irecep]][[jtrait]])<=StatCaseControl[i,j]))/length(lNullRecepDiff[[irecep]][[jtrait]])
	    }
	}
    }
    
    rownames(DiffpvalAbs)<-rownames(StatCaseControl)
    colnames(DiffpvalAbs)<-colnames(StatCaseControl)
    DiffpvalAbsBonf<-DiffpvalAbs*ncol(DiffpvalAbs)*nrow(DiffpvalAbs)
    rownames(DiffpvalRightTail)<-rownames(StatCaseControl)
    colnames(DiffpvalRightTail)<-colnames(StatCaseControl)
    DiffpvalRightTailBonf<-DiffpvalRightTail*ncol(DiffpvalRightTail)*nrow(DiffpvalRightTail)
    rownames(DiffpvalLeftTail)<-rownames(StatCaseControl)
    colnames(DiffpvalLeftTail)<-colnames(StatCaseControl)
    DiffpvalLeftTailBonf<-DiffpvalLeftTail*ncol(DiffpvalLeftTail)*nrow(DiffpvalLeftTail)
    
    save(lNullRecepDiff,StatCaseControl,DiffpvalAbs,DiffpvalAbsBonf,DiffpvalRightTail,DiffpvalRightTailBonf,DiffpvalLeftTail,DiffpvalLeftTailBonf,file=outfile)
}
