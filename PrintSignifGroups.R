PrintSignifGroups<-function(disease,filename,pvalcutoff=0.1){
    #cat(paste(disease,filename))
    load(filename)
    for (lgroup in lsignifGroups){
	if (!is.null(lgroup$pval)){
	    if (lgroup$pval<pvalcutoff){
		#outtext<-""
		for (i in 1:length(lgroup$vars)){
		    for (j in 1:ncol(lgroup$vars[[i]])){
			outtext<-paste(names(lgroup$vars)[i]," ",lgroup$vars[[i]][1,j]," significantly",sep="")
			if (lgroup$vars[[i]][2,j]=="left"){outtext<-paste(outtext," LESSER than Normal.",sep="")}
			if (lgroup$vars[[i]][2,j]=="right"){outtext<-paste(outtext," GREATER than Normal.",sep="")}
			if (lgroup$vars[[i]][2,j]=="abs"){outtext<-paste(outtext," DIFFERENT than Normal.",sep="")}		
			#print(outtext)
			cat(outtext)
			cat("\n")			
		    }
		}	
		#print(lgroup$vars)
		#print(paste("p-value:",lgroup$pval))
		cat(paste("p-value:",lgroup$pval))
		cat("\n")			
		#print("======================================================================================================================")    
		cat("======================================================================================================================")    
		cat("\n")			
	    }
	}
    }
}
