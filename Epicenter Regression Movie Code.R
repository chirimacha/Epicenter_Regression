#===========================================
#   This workspace visualizes the MCMC fitting process of epicenter regression.
# 	It is meant to be run subsequent to the Epicenter_Regression_Model_Simulated_Data.r code.
#===========================================

chains_data<-CHAINS1

#===========================================
#Calculates the number of cases in each age group.
#===========================================
 
max_age<-90

cases_by_age<-rep(0,max_age)
for (w in min(age):max_age)
{
    byAge<-which(age==w)
    cases_by_age[w]<-length(which(POS[byAge]==1))
}
count<-rep(0,max_age)
for (w in min(age):max_age)
{
    count[w]<-length(which(age==w))
}

pcount<-count

#===========================================
#    Function to extract mean estimated prevalence for each year of age
#===========================================

getageprob<-function(a1,a2,model_prob)
{
	cases<-sum(model_prob[age>=a1&age<=a2],na.rm=T)
	count<-length(is.na(model_prob[age>=a1&age<=a2])==F)
	return(list(cases,count))
}


#===========================================
#    Movie function
#    -takes a chain 'ch' from the epicenter regression
#    -The number of epicenters in the model 'epicenters'
#    -'first'  The spot in the chain to begin the movie
#    -'last'  The spot in the chain to stop the movie
#===========================================
#plotting parameters
widthT<-c(.5,.5,.5,.5,.5,.5,.5,.5,.5,0,0)          

MOVIE<-function(ch,epicenters,first,last)
{
	
    reps<-(last-first)+1
	
#set up the layout of the movie
    nf <- layout(matrix(c(1,1,2,2,1,1,3,3,1,1,4,4,1,1,5,5,1,1,6,6), 5, 4, byrow=TRUE),width=c(8,8,4,4),height=c(8,4,4,4,4),respect=TRUE)
	
	epicenters<<-epicenters      
	endT<-epicenters
	endE<-endT+epicenters
	endRF<-endE+2+num_covariates
	
#Array of distance to each epicenter
    DIS<-matrix(nrow = length(lat), ncol = epicenters)
	
#Array of exposure time
    EXP<-matrix(nrow = length(lat), ncol = epicenters)
	
#vector for maximum of exposure time from all epicenter or age (if younger than exposure time)
    maxEXP<-1:length(lat)*0  
    Tchain<-matrix(nrow = reps, ncol = epicenters)
    Echain<-matrix(nrow = reps, ncol = epicenters)
	
	
#Loops through the chain, taking the values for each parameter
	for (u in first:last)
	{
		
#Assign parameter values
		R<-ch$Rchain[u]
		B<-ch$Bchain[u]
		L<-ch$Lchain[u]
        IntroT<-ch[u,1:endT]
		epi<-ch[u,(endT+1):endE]
		RF<-ch[u,(endE+3):endRF]
		epi<-unlist(epi)
		RF<-unlist(RF)
		IntroT<-unlist(IntroT)
		
# Update the exposure matrix and exposure time for each individual
		for ( i in 1:epicenters)
		{
			DIS[,i]<-dise(epi[i])
			EXP[,i]<-(IntroT[i]-(DIS[,i]/R))
		}
		maxEXP<-get_maxEXP(EXP)
		
# Calculate the probability of infection
		model_prob<-getprob(maxEXP,B,RF)
		
#===========================================
#    Create a Map of each observation
#        -Circles are scaled by exposure time
#        -Infected individuals in red
#        -Time of introduction for each micro-epidemic is printed
#===========================================
		
		plot(long,lat,cex=maxEXP+0.5, axes=FALSE)
        axis(1, labels=FALSE)
		axis(2, labels=FALSE)
		box(lty="solid")  
		title(main=paste("Chain", u, "     Likelihood = ",signif(L,digits=3),sep=" "),cex.main=2,outer=FALSE)
		points(long[POS==1],lat[POS==1],cex=maxEXP[POS==1],col="red",pch=21,lwd=4)
		
#Print introduction time of each micro-epidemic
		for ( i in 1:epicenters)
		{
			points(ulong[epi[i]],ulat[epi[i]],cex=9,bg=i+2, pch=23)
			text(ulong[epi[i]],ulat[epi[i]], labels=round(IntroT[i]),cex=1.8)
		}
		
#===========================================
#    Plots Age prevalence curves
#        -Observed in Black
#        -Model estimate in Blue
#===========================================
		
		estimated_cases_by_age<-rep(0,max_age)
		for (i in 1:max_age)
		{
			estimated_cases_by_age[i]<-getageprob(i,i,model_prob)[[1]]
			pcount[i]<-getageprob(i,i,prob)[[2]]
		}
		
#plot real observed data
		plot(cases_by_age/count,cex=count/10,ylim=c(0,.15),col=1,pch=21,xlab="Age",ylab="Prevalence",cex.axis=1.4,lwd=2)
#plot model estimate for chain u, plotted in blue
		points(estimated_cases_by_age/pcount,cex=pcount/10,ylim=c(0,.15),col=4,pch=21,lwd=2)
		
#===========================================
#    Barplots of parameter values
#    -Times of introduction
#    -Rate of spred
#    -Effect of covariates
#===========================================              
		
#Time of introduction
		INTROT<-unlist(IntroT)
		barplot(horiz=TRUE,height=INTROT[1],xlim=c(0,40),col=3,xlab="Time of Introduction",                     cex.axis=2,cex=2,width=widthT[1],ylim=c(0,4))
		
		barplot(horiz=TRUE,height=INTROT[2],xlim=c(0,40),col=4,xlab="Time of Introduction",                      cex.axis=2,cex=2,width=widthT[2],ylim=c(0,4))
		
#Rate of spread
		barplot(horiz=TRUE,height=R,xlim=c(0,80),col="orange",xlab="RATE", cex.axis=2,cex=2,names.arg="  Rate")
        
		RFcol<-c(0,0)
		ifelse(RF[1]>0, RFcol[1]<-"black",RFcol[1]<-"red")
		ifelse(RF[2]>0, RFcol[2]<-"black",RFcol[2]<-"red")
		
#Effect of Risk Factor as a risk ratio
		barplot(horiz=TRUE,height=exp(RF[1]),xlim=c(exp(-.06),exp(.06)),col=RFcol[1],xlab="Increase of Risk", ylab="",xpd=FALSE,cex.axis=1.4,cex=2,names.arg="    Risk Factor")

#allows control of the speed of the movie by requiring a keystroke between frames
		readline("")
	}
}

#=================================
#Open a window for graphics
#Use quartz for Macs
#Use window for Windows
#==================================

quartz(width=10.8,height=8.8)
#windows(width=12, height=12, pointsize=8)

#=================================
#Run the movie function for
#    -chain chains_data
#    -2 epicenters
#=================================

MOVIE(chains_data, 2, 1, 100)

#======================================================
#Determine how often each house is chosen as an epicenter
#	-compare to the true values from the simulation (house 100 and house 780)
#======================================================

epiFreq<-rep(0,length(lat))
E1<-chains_data$E1
E2<-chains_data$E2
for (h in 1:length(lat))
{
    epiFreq[h]<-length(which(E1==h))+length(which(E2==h))
}
par(mfrow=c(1,2))
plot(long,lat, cex=epiFreq+0.2, main="Posterior Estimate of Location of Epicenters")
plot(long,lat, cex=0.2, main="True location of Epicenters")

#plot the location of the true epicenters
example_epi<-c(100,780)
points(long[example_epi],lat[example_epi],cex=5,col=3, pch=19)

#=====================================
#Examine the chains using the coda package
#===================================== 

#library(coda)
#?coda
