library(mvSLOUCH)
library(permute)
source("wholePipeline.R")
source("ModelAnalysis.R") 
source("SimulateHealthyModel.R")
source("getNullDistribution.R")
source("CompareToNull.R")
source("CompareToNullGroups.R")
source("FindSignif.R")
source("PrintSignifGroups.R")


disease<-"Copd" ## prefix of disease, the data files have to be called diseaseNormal.csv and diseasePatient.csv (like you have now)

treefile<-"receptorsMrBayes.tre"
varnames<-c("Exp","Meth") ## the names of the columns in the data files corresponding to our variables of interest, 
## can be any number so you can use this for cancer and diabetes. The names have to be the same for the Normal and Patients data files!

Merrorvarnames<-c("varExp","varMeth") ## the names of the columns in the data files corresponding to the measurement errors of our  variables of interest, 
## this has to be the same length as varnames! The order also has to be the same! The code assumes that the name of the i-th measurement error column
## in Merrorvarnames corresponds to the i-th variable in varnames. Therefore if you will have a variable without measurement error you need
## to insert a dummy variable into your data file.
## The names have to be the same for both Normal and Patient data files!

simset<-1
numSims<-200000
pvalcutoff<-0.05
maxnumsigs<-6 ## max number of variables in group

directCompare<-FALSE ## this is just which way of assessing significance we want
IndepDiff<-FALSE  ## this is just which way of assessing significance we want
## The below TRUE/FALSE vector tells us which parts of the analysis to run, as all take a long time we may just run more downstream ones
## using previously done computations. Of course if we have a FALSE between TRUE values then the results of the first TRUE will not 
## have any effect on the results of the second TRUE. Each stage of the analysis is based solely on the results of the previous stage.
## The elements of the vtorun vector are the following :
## vtorun[1] : should we estimate the model parameters and find the best model
## vtorun[2] : should we simulate data under the best found model for the Normal group
## vtorun[3] : should we extract the desired (direct comparison, conditional on healthy, independent differences) null distribution from the simulated data
## vtorun[4] : should we compare our true data to the extracted null distribution
## vtorun[5] : should we extract the significantly different (lying in either tail) receptors from the comparison between true data and null distribution
## vtorun[6] : should we extract the significantly different groups (composed of individually significantly different) of receptors from the comparison between true data and null distribution
vtorun<-c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
##now run the analysis

mvslPipeline(disease=disease,treefile=treefile,varnames=varnames,Merrorvarnames=Merrorvarnames,vtorun=vtorun,numSims=numSims,directCompare=directCompare,IndepDiff=IndepDiff,simset=simset,pvalcutoff=pvalcutoff,maxnumsigs=maxnumsigs)

directCompare<-TRUE  ## this is just which way of assessing significance we want
IndepDiff<-FALSE  ## this is just which way of assessing significance we want
## The below TRUE/FALSE vector tells us which parts of the analysis to run, as all take a long time we may just run more downstream ones
## using previously done computations. Of course if we have a FALSE between TRUE values then the results of the first TRUE will not 
## have any effect on the results of the second TRUE. Each stage of the analysis is based solely on the results of the previous stage.
## The elements of the vtorun vector are the following :
## vtorun[1] : should we estimate the model parameters and find the best model
## vtorun[2] : should we simulate data under the best found model for the Normal group
## vtorun[3] : should we extract the desired (direct comparison, conditional on healthy, independent differences) null distribution from the simulated data
## vtorun[4] : should we compare our true data to the extracted null distribution
## vtorun[5] : should we extract the significantly different (lying in either tail) receptors from the comparison between true data and null distribution
## vtorun[6] : should we extract the significantly different groups (composed of individually significantly different) of receptors from the comparison between true data and null distribution
## Notice that here we have no need to estimate the best model as it was already done, here we just use a different definition of the null
vtorun<-c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE)
##now run the analysis

mvslPipeline(disease=disease,treefile=treefile,varnames=varnames,Merrorvarnames=Merrorvarnames,vtorun=vtorun,numSims=numSims,directCompare=directCompare,IndepDiff=IndepDiff,simset=simset,pvalcutoff=pvalcutoff,maxnumsigs=maxnumsigs)

directCompare<-FALSE  ## this is just which way of assessing significance we want
IndepDiff<-TRUE  ## this is just which way of assessing significance we want
vtorun<-c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE)
## The below TRUE/FALSE vector tells us which parts of the analysis to run, as all take a long time we may just run more downstream ones
## using previously done computations. Of course if we have a FALSE between TRUE values then the results of the first TRUE will not 
## have any effect on the results of the second TRUE. Each stage of the analysis is based solely on the results of the previous stage.
## The elements of the vtorun vector are the following :
## vtorun[1] : should we estimate the model parameters and find the best model
## vtorun[2] : should we simulate data under the best found model for the Normal group
## vtorun[3] : should we extract the desired (direct comparison, conditional on healthy, independent differences) null distribution from the simulated data
## vtorun[4] : should we compare our true data to the extracted null distribution
## vtorun[5] : should we extract the significantly different (lying in either tail) receptors from the comparison between true data and null distribution
## vtorun[6] : should we extract the significantly different groups (composed of individually significantly different) of receptors from the comparison between true data and null distribution
## Notice that here we have no need to estimate the best model as it was already done, here we just use a different definition of the null
##now run the analysis
mvslPipeline(disease=disease,treefile=treefile,varnames=varnames,Merrorvarnames=Merrorvarnames,vtorun=vtorun,numSims=numSims,directCompare=directCompare,IndepDiff=IndepDiff,simset=simset,pvalcutoff=pvalcutoff,maxnumsigs=maxnumsigs)

