FindSignif<-function(disease,filename,pvalcutoff=0.1){
    cat(paste(disease,filename))
    cat("\n");    cat("\n")
    load(filename)
    varnames<-colnames(DiffpvalAbs)
    vpvalsabs<-sort(c(DiffpvalAbs),na.last=FALSE)
    vpvalslefttail<-sort(c(DiffpvalLeftTail),na.last=FALSE)
    vpvalsrighttail<-sort(c(DiffpvalRightTail),na.last=FALSE)
    #lsignifabs<-vector("list",length(varnames));names(lsignifabs)<-varnames
    #lsigniflefttail<-vector("list",length(varnames));names(lsigniflefttail)<-varnames
    #lsignifrighttail<-vector("list",length(varnames));names(lsignifrighttail)<-varnames
    lsignifsFound<-vector("list",length(varnames));names(lsignifsFound)<-varnames
    for (v in varnames){
	dpvalabs<-which(DiffpvalAbs[,v]<pvalcutoff)
	dpvalabs<-setdiff(dpvalabs,which(is.na(DiffpvalAbs[,v])))
	dpvalleft<-which(DiffpvalLeftTail[,v]<pvalcutoff)
	dpvalleft<-setdiff(dpvalleft,which(is.na(DiffpvalLeftTail[,v])))
	dpvalright<-which(DiffpvalRightTail[,v]<pvalcutoff)
	dpvalright<-setdiff(dpvalright,which(is.na(DiffpvalRightTail[,v])))
	    if (length(dpvalabs)>0){
		cat(c(disease,v,"Significantly DIFFERENT than Normal, p-values"))
		cat("\n")
		res<-matrix(DiffpvalAbs[dpvalabs,v],nrow=1)
		rownames(res)<-v
		colnames(res)<-rownames(DiffpvalAbs)[dpvalabs]
		cat(colnames(res))
		cat("\n")
		cat(res)
		cat("\n")
		#lsignifabs$v<-colnames(res)
		i<-which(names(lsignifsFound)==v)
		lsignifsFound[[i]]<-rbind(colnames(res),"abs")
		rownames(lsignifsFound[[i]])<-c("receptors","tail")
		cat("======================================================================================================================")
	        cat("\n")
	    }
	    NAnum<-length(which(is.na(DiffpvalAbs)))
	    m<-length(vpvalsabs)-NAnum
	    vBHvals<-pvalcutoff*(1:m)/m
	    vBHp<-c()
	    if (NAnum==0){vBHp<-which(vpvalsabs-vBHvals <=0)}
	    else{vBHp<-which(vpvalsabs[-which(is.na(DiffpvalAbs))]-vBHvals <=0)}
	    if (length(vBHp)>0){
		pvalcutoffBH<-vpvalsabs[max(vBHp)]
	    	dpvalabs<-which(DiffpvalAbs[,v]<pvalcutoffBH)
	    	dpvalabs<-setdiff(dpvalabs,which(is.na(DiffpvalAbs[,v])))
	    	if (length(dpvalabs)>0){
	    	    cat(c(disease,v,"Significantly DIFFERENT than Normal, BH controlled p-values"))
		    cat("\n")
		    res<-matrix(DiffpvalAbs[dpvalabs,v],nrow=1)
		    rownames(res)<-v
		    colnames(res)<-rownames(DiffpvalAbs)[dpvalabs]
		    cat(colnames(res))
		    cat("\n")
		    cat(res)
		    cat("\n")
		    cat("======================================================================================================================")
		    cat("\n")
		    cat("======================================================================================================================")
		    cat("\n")
		}
	    }
	    if (length(dpvalleft)>0){
		cat(c(disease,v,"Significantly LESSER than Normal, p-values"))
		cat("\n")
		res<-matrix(DiffpvalLeftTail[dpvalleft,v],nrow=1)
		rownames(res)<-v
		colnames(res)<-rownames(DiffpvalLeftTail)[dpvalleft]
		cat(colnames(res))
		cat("\n")
		cat(res)
		cat("\n")
		i<-which(names(lsignifsFound)==v)
		if (is.null(lsignifsFound[[i]])){
		    lsignifsFound[[i]]<-rbind(colnames(res),"left")
		    rownames(lsignifsFound[[i]])<-c("receptors","tail")
		}else{lsignifsFound[[i]]<-cbind(lsignifsFound[[i]],rbind(colnames(res),"left"))}
		#lsigniflefttail$v<-colnames(res)
		cat("======================================================================================================================")
	        cat("\n")
	    }
	    NAnum<-length(which(is.na(DiffpvalLeftTail)))
	    m<-length(vpvalslefttail)-NAnum
	    vBHvals<-pvalcutoff*(1:m)/m    
	    vBHp<-c()
	    if (NAnum==0){vBHp<-which(vpvalslefttail-vBHvals <=0)}
	    else{vBHp<-which(vpvalslefttail[-which(is.na(DiffpvalLeftTail))]-vBHvals <=0)}
	    if (length(vBHp)>0){
		pvalcutoffBH<-vpvalslefttail[max(vBHp)]
	    	dpvalleft<-which(DiffpvalLeftTail[,v]<pvalcutoffBH)
	    	dpvalleft<-setdiff(dpvalleft,DiffpvalLeftTail[,v])
	    	if (length(dpvalleft)>0){
	    	    cat(c(disease,v,"Significantly LESSER than Normal, BH controlled p-values"))
		    cat("\n")
		    res<-matrix(DiffpvalLeftTail[dpvalleft,v],nrow=1)
		    rownames(res)<-v
		    colnames(res)<-rownames(DiffpvalLeftTail)[dpvalleft]
		    cat(colnames(res))
		    cat("\n")
		    cat(res)
		    cat("\n")
		    cat("======================================================================================================================")
		    cat("\n")
		    cat("======================================================================================================================")
		    cat("\n")
		}
	    }	
	    if (length(dpvalright)>0){
		cat(c(disease,v,"Significantly GREATER than Normal, p-values"))
		cat("\n")
		res<-matrix(DiffpvalRightTail[dpvalright,v],nrow=1)
		rownames(res)<-v
		colnames(res)<-rownames(DiffpvalRightTail)[dpvalright]
		cat(colnames(res))
		cat("\n")
		cat(res)
		cat("\n")
		i<-which(names(lsignifsFound)==v)
		if (is.null(lsignifsFound[[i]])){
		    lsignifsFound[[i]]<-rbind(colnames(res),"right")
		    rownames(lsignifsFound[[i]])<-c("receptors","tail")
		}else{lsignifsFound[[i]]<-cbind(lsignifsFound[[i]],rbind(colnames(res),"right"))}
		#lsignifrighttail$v<-colnames(res)
		cat("======================================================================================================================")
	        cat("\n")
	    }
	    NAnum<-length(which(is.na(DiffpvalRightTail)))
	    m<-length(vpvalsrighttail)-NAnum
	    vBHvals<-pvalcutoff*(1:m)/m    
	    vBHp<-c()
	    if (NAnum==0){vBHp<-which(vpvalsrighttail-vBHvals <=0)}
	    else{vBHp<-which(vpvalsrighttail[-which(is.na(DiffpvalRightTail))]-vBHvals <=0)}
	    if (length(vBHp)>0){
		pvalcutoffBH<-vpvalsrighttail[max(vBHp)]
	    	dpvalright<-which(DiffpvalRightTail[,v]<pvalcutoffBH)
	    	dpvalright<-setdiff(dpvalright,which(is.na(DiffpvalRightTail[,v])))
		if (length(dpvalright)>0){
		    cat(c(disease,v,"Significantly GREATER than Normal, BH controlled p-values"))
		    cat("\n")
		    res<-matrix(DiffpvalRightTail[dpvalright,v],nrow=1)
		    rownames(res)<-v
		    colnames(res)<-rownames(DiffpvalRightTail)[dpvalright]
	    	    cat(colnames(res))
	    	    cat("\n")
	    	    cat(res)
		    cat("\n")
		    cat("======================================================================================================================")
		    cat("\n")
		    cat("======================================================================================================================")
		    cat("\n")
		}
	    }
	cat("**********************************************************************************************************************")
	cat("\n")
    }
    list("SignificantFound"=lsignifsFound)#"Abs"=lsignifabs,"LeftTail"=lsigniflefttail,"RightTail"=lsignifrighttail)
}
