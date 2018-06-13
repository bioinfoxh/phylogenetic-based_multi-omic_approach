source("EvolModelEst.R")

FindBestModel<-function(treefile,datafile,varnames,Merrorvarnames,outprefix,numreps=3,regimes=NULL){
## varnames and Merrorvarnames have to be in the same order!!
## regimes is not needed at the moment it kept for future development!

    ## preparing the data for the mvSLOUCH package
    TraitData<-read.csv(datafile,header=TRUE,row.names=1)
    apetree<-read.nexus(treefile)
    apetree<-drop.tip(apetree,setdiff(apetree$tip.label,rownames(TraitData)))
    ouchtree<-ape2ouch(apetree,scale=1)
    ouchtree@nodelabels[1:(ouchtree@nnodes-ouchtree@nterm)]<-as.character(1:(ouchtree@nnodes-ouchtree@nterm))

    TraitData<-TraitData[as.character(ouchtree@nodelabels[(ouchtree@nnodes-ouchtree@nterm+1):ouchtree@nnodes]),]

##    regimesFitch<-fitch.mvsl(ouchtree,regimes,root=1,deltran=TRUE) ## we do not need this at the moment but kept just in case for future development

    dfdata<-TraitData[,varnames]
    numtraits<-length(varnames)
## --------------------------------------------------------------

    M.error<-rep(0,ncol(dfdata)*nrow(TraitData))
     if (!is.null(Merrorvarnames)){
         for (i in 1:length(Merrorvarnames)){
    		M.error[which((1:length(M.error))%%numtraits==(i%%length(Merrorvarnames)))]<-TraitData[,Merrorvarnames[i]]
        }
    }
    M.error<-diag(M.error)

    lBestModel<-evol.model.est(ouchtree,dfdata,regimes=NULL,root.regime=NULL,M.error=M.error,repeats=3,model.setups=NULL,predictors=2,kY=1,doPrint=TRUE)
    BestModel<-lBestModel$BestModel
    sink(paste(outprefix,"_BestModel.txt",sep=""))
    print(lBestModel$BestModel)
    sink()

    save(lBestModel,BestModel,file=paste(outprefix,"_BestModel.RData",sep=""))
    save.image(file=paste(outprefix,"_mvslouch.RData",sep=""))
}
