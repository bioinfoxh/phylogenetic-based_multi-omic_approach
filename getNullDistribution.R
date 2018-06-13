getNullDistribution<-function(controldatafile,casedatafile,treefile,simulateddatafile,varnames,Merrorvarnames=NULL,outfile,directCompare=FALSE,IndepDiff=FALSE){
## order of varnames and Merrorvarnames has to be the same

    load(simulateddatafile) ## the variable has to have a specific name!!
    numtraits<-ncol(lSimulatedHealthyData[[1]])

    Controldata<-read.csv(controldatafile,header=TRUE,row.names=1)
    Receptorstreeape<-read.nexus(treefile)
##    Receptorstreeape<-drop.tip(Receptorstreeape,setdiff(Receptorstreeape$tip.label,rownames(Controldata)))
    Receptorstreeouch<-ape2ouch(Receptorstreeape,scale=1)
    Receptorstreeouch@nodelabels[1:(Receptorstreeouch@nnodes-Receptorstreeouch@nterm)]<-as.character(1:(Receptorstreeouch@nnodes-Receptorstreeouch@nterm))
    numtips<-phyltreeouch@nterm
    numinternal<-Receptorstreeouch@nnodes-Receptorstreeouch@nterm
    ouchtree<-Receptorstreeouch
    Controldata<-Controldata[as.character(Receptorstreeouch@nodelabels[(Receptorstreeouch@nnodes-Receptorstreeouch@nterm+1):Receptorstreeouch@nnodes]),]
    BaseData<-Controldata


    dfControldata<-Controldata[,varnames]
    if (length(varnames)==1){
	dfControldata<-matrix(dfControldata,ncol=1,nrow=length(dfControldata))
        colnames(dfControldata)<-varnames
    }

    M.error<-rep(0,ncol(dfControldata)*nrow(Controldata))
    if (!is.null(Merrorvarnames)){
        for (i in 1:length(Merrorvarnames)){
	    M.error[which((1:length(M.error))%%numtraits==(i%%length(Merrorvarnames)))]<-Controldata[,Merrorvarnames[i]]
	}
    }
    M.error[which(is.na(M.error))]<-0
    varControl<-M.error

    M.error<-NA
    Casedata<-read.csv(casedatafile,header=TRUE,row.names=1)
    Casedata<-Casedata[as.character(Receptorstreeouch@nodelabels[(Receptorstreeouch@nnodes-Receptorstreeouch@nterm+1):Receptorstreeouch@nnodes]),]
    dfCasedata<-Casedata[,varnames]
    if (length(varnames)==1){
        dfCasedata<-matrix(dfCasedata,ncol=1,nrow=length(dfCasedata))
	colnames(dfCasedata)<-varnames
    }
    M.error<-rep(0,ncol(dfCasedata)*nrow(Casedata))
    if (!is.null(Merrorvarnames)){
        for (i in 1:length(Merrorvarnames)){
		M.error[which((1:length(M.error))%%numtraits==(i%%length(Merrorvarnames)))]<-Casedata[,Merrorvarnames[i]]
	}
    }
    M.error[which(is.na(M.error))]<-0
    varCase<-M.error

    BaseData<-BaseData[,varnames]
    if (directCompare){BaseData<-matrix(0,ncol=ncol(BaseData),nrow=nrow(BaseData))}
    if (length(varnames)==1){BaseData<-matrix(BaseData,ncol=length(varnames),nrow=length(BaseData))}
    colnames(BaseData)<-varnames

    if (!IndepDiff){
        nullvalues<-sapply(lSimulatedHealthyData,function(x,y,numinternal,varsused,varCase)
	{((x[-(1:numinternal),varsused]
	+matrix(rmvnorm(1,mean=rep(0,length(varCase)),sigma=diag(varCase)),ncol=length(varsused),nrow=nrow(x)-numinternal,byrow=TRUE))
	-(y))},
	y=BaseData,numinternal=numinternal,varsused=varnames,varCase=varCase,simplify=FALSE)
    }else{
	nullvalues<-sapply(1:(length(lSimulatedHealthyData)%/%2),function(i,lData,numinternal,varCase,varControl){
	x<-lData[[2*(i-1)+1]]; y<-lData[[2*i]][-(1:numinternal),];
	((x[-(1:numinternal),]+matrix(rmvnorm(1,mean=rep(0,length(varCase)),sigma=diag(varCase)),ncol=ncol(x),nrow=nrow(x)-numinternal,byrow=TRUE))
	-(y+matrix(rmvnorm(1,mean=rep(0,length(varControl)),sigma=diag(varControl)),ncol=ncol(y),nrow=nrow(y),byrow=TRUE))
	)},lData=lSimulatedHealthyData,numinternal=numinternal,varCase=varCase,varControl=varControl,simplify=FALSE)
    }
    
    lNullRecepDiff<-vector("list",numtips)
    names(lNullRecepDiff)<-rownames(lSimulatedHealthyData[[1]])[-(1:numinternal)]
    lNullRecepDiff<-sapply(lNullRecepDiff,function(x,k,vnames){
		ltmp<-vector("list",k);names(ltmp)<-vnames;for (i in 1:numtraits){ltmp[[i]]<-c(NA)};ltmp}
		,k=numtraits,vnames=colnames(lSimulatedHealthyData[[1]]),simplify=FALSE)

    for (i in 1:length(nullvalues)){
	for (j in 1:numtips){
	    for (k in 1:numtraits){
		lNullRecepDiff[[j]][[k]]<-c(lNullRecepDiff[[j]][[k]],nullvalues[[i]][j,k])
	    }
	}
    }
    lunsortedNullRecepDiff<-lNullRecepDiff
    for (j in 1:numtips){
	for (k in 1:numtraits){
	    lNullRecepDiff[[j]][[k]]<-sort(lNullRecepDiff[[j]][[k]][-1],na.last=FALSE)
	}	
    }
    save(lNullRecepDiff,lunsortedNullRecepDiff,nullvalues,file=outfile)
}
