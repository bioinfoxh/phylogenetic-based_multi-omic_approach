CompareToNullGroups<-function(controldatafile,casedatafile,nulldistributionfile,varnames,lUniVars,outfile,directCompare=FALSE,maxnumsigs=2){
    Controldata<-read.csv(controldatafile,header=TRUE,row.names=1)
    Casedata<-read.csv(casedatafile,header=TRUE,row.names=1)
    load(nulldistributionfile) ## here we force the variable name to be lNullRecepDiff
    if (directCompare){StatCaseControl<-Casedata[,varnames]}
    else{StatCaseControl<-(Casedata[,varnames]-Controldata[,varnames])}
    if (length(varnames)==1){
	StatCaseControl<-matrix(StatCaseControl,ncol=1,nrow=length(StatCaseControl))
	colnames(StatCaseControl)<-varnames
    }
    
    lcombNums<-vector("list",length(lUniVars))
    k<-1
    for (i in 1:length(lUniVars)){
	if (!is.null(lUniVars[[i]])){ ## if no receptor was significant NULL returned
	    for(j in 1:ncol(lUniVars[[i]])){
		if (lUniVars[[i]][2,j]!="abs"){ ## ignore the abs ones we only want up or down regulated
		    lcombNums[[k]]<-c(i,j); k<-k+1
		}	
    	    }
	}
    }
    

    numsigs<-length(lcombNums)
    lcombsindex<-list()
    k<-1
    for (i in 2:min(numsigs,maxnumsigs)){    
	lcombsindex[[k]]<-combn(1:numsigs,i)
	k<-k+1
    }
    lcombsAlltmp<- sapply(lcombsindex,function(x,lcombNums){
	    res<-sapply( split(x, rep(1:ncol(x), each = nrow(x))),function(y,lcombNums){
		res<-NULL
		bNULL<-FALSE
		for (z in y){if (is.null(lcombNums[[z]])){bNULL<-TRUE}}
		if (!bNULL){res<-sapply(lcombNums[y],function(z){z})}
		res
	    },lcombNums=lcombNums,simplify=FALSE)
	    names(res)<-NULL
	    res
	},lcombNums=lcombNums,simplify=FALSE )

lcombsAlltmp.2<-list()
for (lcomb in lcombsAlltmp){lcombsAlltmp.2<-c(lcombsAlltmp.2,lcomb)}
    lcombsAll<- sapply(lcombsAlltmp.2,function(x,lUniVars){
	res<-NA
	if (!is.null(x)){
	    signifvars<-unique(x[1,])
	    res<-vector("list",length(signifvars))
	    names(res)<-names(lUniVars)[signifvars]
	    for (v in names(res)){
		vindex<-which(names(lUniVars)==v)
		r<-which(names(res)==v)
		res[[r]]<-matrix(NA,nrow=2,ncol=length(which(x[1,]==vindex)))
		for (j in 1:ncol(res[[r]])){
		    k<-x[2,which(x[1,]==vindex)[j]]
    		    ru<-which(names(lUniVars)==v)
		    res[[r]][,j]<-lUniVars[[ru]][,k]
		}	    
	    }
	}
	res
    },lUniVars=lUniVars,simplify=FALSE)

    ##generate from lUniVars all possible combinations
#we have 3 list abs, left, right
#abs should be ignored only left, right
#use lunsortedNullRecepDiff and then take the intersection of tests first element is NA!!
#check if lunsorted really has the order kept

    lsignifGroups<-vector("list",length(lcombsAll))
    k<-1
    for (lcomb in lcombsAll){
	if ((!is.null(lcomb))&&(!is.na(lcomb[1]))){
	    vsignifcombs<-1:length(lunsortedNullRecepDiff[[1]][[1]]) ## should be OK we have same sample length eveyrwhere

    #print("this needs to be checked how it is done")
	    for (i in 1:length(lcomb)){    
    		for (j in 1:ncol(lcomb[[i]])){## first row is receptor second is abs/right/left
    	    	    irecep<-which(names(lunsortedNullRecepDiff)==lcomb[[i]][1,j]) ##rownames(StatCaseControl)[i])
		    jtrait<-which(names(lunsortedNullRecepDiff[[irecep]])==names(lcomb)[i]) ##colnames(StatCaseControl)[j])
		    vcurrsigs<-c()
		    iSCC<-which(rownames(StatCaseControl)==names(lunsortedNullRecepDiff)[irecep])
		    jSCC<-which(colnames(StatCaseControl)==names(lunsortedNullRecepDiff[[irecep]])[jtrait])
		    if ((!is.na(StatCaseControl[iSCC,jSCC]))&&(!is.na(Casedata[iSCC,jSCC]))&&(!is.na(Controldata[iSCC,jSCC]))){
			## here we are doing jointly for many variables so we first calculate the intersection of those null that are less extreme
			## than the observation and then to get the p-value we take 1 - this
			## otherwise we would need to do an alternating sum
			if (lcomb[[i]][2,j]=="abs"){vcurrsigs<-which(abs(lunsortedNullRecepDiff[[irecep]][[jtrait]])<=StatCaseControl[iSCC,jSCC])}#)/length(lunsortedNullRecepDiff[[irecep]][[jtrait]])}
			if (lcomb[[i]][2,j]=="right"){vcurrsigs<-which((lunsortedNullRecepDiff[[irecep]][[jtrait]])<=StatCaseControl[iSCC,jSCC])}#/length(lunsortedNullRecepDiff[[irecep]][[jtrait]])}
			if (lcomb[[i]][2,j]=="left"){vcurrsigs<-which((lunsortedNullRecepDiff[[irecep]][[jtrait]])>=StatCaseControl[iSCC,jSCC])}#/length(lunsortedNullRecepDiff[[irecep]][[jtrait]])}
			vsignifcombs<-intersect(vsignifcombs,vcurrsigs)
#			print(length(vsignifcombs))
		    }
		}		
	    }
#	    print("###############")
	    lsignifGroups[[k]]<-vector("list",2)
	    names(lsignifGroups[[k]])<-c("vars","pval")
	    lsignifGroups[[k]]$vars<-lcomb
	    lsignifGroups[[k]]$pval<-1-length(vsignifcombs)/length(lunsortedNullRecepDiff[[1]][[1]]) ## We take 1- to get the p-value  !!
	    k<-k+1
	}
    }
    save(lunsortedNullRecepDiff,StatCaseControl,lsignifGroups,file=outfile)
}
