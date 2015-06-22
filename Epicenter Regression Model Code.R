#===============================================#
#    Load Data
#    -Data file with covariates and outcome variable
#    -Distance matrix between households, hDIST (potential sites of introduction of the disease agent)
#    -Distance matrix between individuals, pDIST (units of observation) and potential sites of introduction
#===============================================#
simdata<-read.csv("Simulated_data.csv")
attach(simdata)
ulat<-lat
ulong<-long
distance_grid<-as.matrix(dist(cbind(lat,long)))
hDIST<-distance_grid
pDIST<-distance_grid


#===============================================#
#    General Parameters describing the data
#===============================================#
num_obs<-length(lat)
num_houses<-length(lat)
num_covariates<-1
epicenters<-2

#===============================================#
#    General Parameters for the Monte Carlo Markov Chains
#            -length of chains
#            -Thinning
#            -Number of Covariates
#            -Initial Values     
#===============================================#

reps<-1000
thin<-10
chain_length<-reps/thin
THIN<-1:chain_length*10
it<-1
trunc<-1
#Initial Values for the Parameters
init_B=-4
init_R=20
ifelse(num_covariates==1,init_RF<-0, init_RF<-array(0,num_covariates))
#Have to 'unlock' t so that we can use it above,
#otherwise R doesnt allow changes to it
t<-0
rf<-0
tnew_lik<-0

#===============================================#
#    Key Model Functions
#		1. Distance of each person's household to an epicenter
#       2. Exposure time
#       3. Probability of infection per unit of time exposed
#       4. Likelihood Function
#===============================================#


#Function to calculate the distance between a person's household and an epicenter
 
 dise<-function(epi)
	{
	return(pDIST[,epi])
	}

# Function that calculates the time that each individual is exposed as the time that her household is exposed
# or her age if she was born into an exposed household.

 get_maxEXP<-function(expt)
    {
     maxexpt<-apply(expt,1,max)
     ind<-1:length(maxexpt)*(age<=maxexpt)
     ind[is.na(ind)]<-0
     maxexpt[ind]<-age[ind]
     maxexpt[maxexpt<0]<-0  #set exposure time to 0 for those the epidemic has not yet reached
     return(maxexpt)
    }

 #Functions to calculate the cumulative probability of infection over time exposed
 #    -Version 1: for a single covariate (RF)
 #     -Version 2: for multiple covariates (RF is a vector)

#Version 1- a single covariate

 if (num_covariates==1)
    {
     getprob<-function(maxEXP,B,RF)
		{
		prob<-1-exp(exp(B+(RF*X))*-maxEXP)
		return(prob)
		}
    }

#Version 2- multiple covariates

 if(num_covariates>1)
    {
     getprob<-function(maxEXP,B,RF)
		{
		covars<-array(0,c(num_obs,num_covariates))
		for(i in 1:num_covariates)
			{
			covars[,i]<-X[[i]]*RF[i]
			}
     covars<-rowSums(covars)
     prob<-1-exp(exp(B+covars)*-maxEXP)
     return(prob)
		}
    }

    #Likelihood function
    #  Calculates the Binomial Likelihood, the constant is omitted as it cancels out when comparing        
	#  likelihoods in the Metropolis algorithm.   

lik<-function(prob,maxexp)
     {
     L<-prod(prob[POS==1])*prod(1-prob[POS==0])
     return(L)
     }

   


#===============================================#
#    MCMC Initialization Functions
#       1. Assignment of Initial Values to parameters
#       2. Creation of empty chains of correct dimensions to store the estimates
#       3. Assignment of parameters of proposal distributions, and counters for acceptance for each parameter
#       4. Get initial Likelihood
#===============================================#

 init_values<-function(epicenters,init_B,init_RF,init_R)

    {
    IntroT<<-runif(epicenters,min=1,max=40)
    B<<-init_B              
	RF<<-init_RF         
    R<<-init_R         
    epi<<-sample(1:num_houses, epicenters, replace=FALSE)
    }

 chain_init<-function(reps,epicenters)
    {
 #Set chains for scalars
		Rchain<<-rep(NA,reps)
		Lchain<<-rep(NA,reps)
		Bchain<<-rep(NA,reps)
 
 #Make a matrix to store the chains of dimensions: # of columns = the number of epicenters, the number of rows=the length of the chain (number of reps)        
		Tchain<<-matrix(nrow = reps, ncol = epicenters)
		Echain<<-matrix(nrow = reps, ncol = epicenters)
 
 #Sets the Risk factor chain as a vector (for 1 risk factor) or a matrix for >1 risk factor
		if(num_covariates==1) {RFchain<<-rep(NA,reps)}
		if(num_covariates>1) {RFchain<<-matrix(NA,nrow = reps, ncol = num_covariates)}
 
 #Array of distance to each epicenter
		DIS<<-matrix(nrow = num_obs, ncol = epicenters)
 
 #Array of exposure time
		EXP<<-matrix(nrow = num_obs, ncol = epicenters)
 
 #Vector for maximum of exposure time from all epicenter or age (if younger than exposure time)
 		maxEXP<<-1:994*0    
	}



proposal_init<-function(epicenters)
    {
 #Variance on proposals for R, the speed of spread of the parasite
		varR<<-80
 #Variance on proposals for B, the risk of infection in the absence of covariates
		varB<<-.25  
 #Variance on proposals for RF, the risk factor(s) for infection given exposure
		if(num_covariates==1) {varRF<<-.001}
		if(num_covariates>1) {varRF<<-c(.001,.1)}
 #Variance on proposals for IntroT, the time(s) of introduction of the parasite into each epicenter
		varT<<-rep(10, epicenters)
 #Variance on proposals for EPI, the location(s) of each epicenter
		varEPI<<-rep(100^2, epicenters)
    }

 
lik_init<-function(epicenters)

    {  
#Initial DIS (distance to epicenter) and EXP (exposure time) matrices
       for ( i in 1:epicenters )
         {
            DIS[,i]<<-dise(epi[i])
            EXP[,i]<<-(IntroT[i]-(DIS[,i]/R))
          }

       maxEXP<<-get_maxEXP(EXP)
       prob<<-getprob(maxEXP,B,RF)
       LIK<<-lik(prob,maxEXP)
	}

#===============================================#
#    Definition of Priors and Functions to return Prior Probability
#    1. Hyperparameters
#       2. Prior mean of R     --Rate of Spread
#       3. Prior mean of B    --Baseline instantaneous Risk of Infection for those without risk factors
#                     mean of this prior corresponds to .5% prevalence among individuals with no insects in
#                      their houses nor animals sleeping inside at night. Variance is very large, making it less informative.
#       4. Prior RF      --Coefficient of the risk factor describing its multiplicative effect on Baseline Risk
#            -the mean is set at 0 (no effect of the risk factor); variance is very large to minimize the  
#             importance of the prior
#       5. Prior IntroT    --Time of Introduction for each epicenter  (Uniform prior)
#       6. Prior EPI     --Location (indicated by house code) of each (Uniform prior)
#===============================================#

prior_mean_R<-20
prior_sd_R<-100^(.5)
prior_meanB<- -5.295812        
prior_sd_B<-10^3                    
prior_meanRF<-0           
prior_sd_RF<-10^3
min_T<-1
max_T<-40

 
R_prior_prob<-function(r)
    {
     prior_prob<-dnorm(r,prior_mean_R,prior_sd_R)
	 return(prior_prob)
    }

B_prior_prob<-function(b)
	{                 
	prior_prob<-dnorm(b,mean=prior_meanB,sd=prior_sd_B)  
	return(prior_prob)
	}

  

RF_prior_prob<-function(rf)
	{
	prior_prob<-dnorm(rf,mean=prior_meanRF,sd=prior_sd_RF)  
	return(prior_prob)	
	}

T_prior_prob<-function(t)
	{
	prior_prob<-dunif(t,min_T,max_T)	
	return(prior_prob)
	}   

EPI_prior_prob<-function(epi)
	{
	prior_prob<-dunif(epi,1,num_houses)
	return(prior_prob)
	}

   

#===============================================#

#    Proposal Probabilities for the Metropolis-Hastings algorithm

#       1. R_proposal     --Normal, centered at the current value, with variance varR

#       2. B_proposal --Normal, centered at the current value with variance varB

#       3. RF_proposal  --Normal, centered at the current value with variance RF_var

#       4. T_proposal     --Normal centered at the current value, truncated if a value <1 or >40 is returned

#       5. EPI_proposal --Normal, centered spatially on the current location.

#===============================================#

R_proposal<-function(r,R)

	{

	proposal_probability<-dnorm(r,R,sd=varR^.5)

    return(proposal_probability)

    }

B_proposal<-function(b,B)

	{

	proposal_probability<-dnorm(b,B,sd=varB^.5)

    return(proposal_probability)

    }

RF_proposal<-function(rf,RF,var_RF)

	{

	sd_RF<-var_RF^.5

	proposal_probability<-dnorm(rf,RF,sd=sd_RF)

	return(proposal_probability)

    }

T_proposal<-function(t,IntroT,var_T)

	{

	sd_T<-var_T^.5

	proposal_probability<-dnorm(t,IntroT,sd=sd_T)

	return(proposal_probability)

    }

EPI_proposal<-function(ep1,ep2,var_EPI)

	{

	sd_EPI<-var_EPI^.5

	weighted_distance<-dnorm(hDIST[,ep1],0,sd=sd_EPI)

	all_distance<-sum(weighted_distance)

	prob_move_ep1_to_ep2<-weighted_distance[ep2]/all_distance

	proposal_probability<-prob_move_ep1_to_ep2

	return(proposal_probability)

    }

#===============================================#
#    MCMC Core Update Functions
#		Each parameter is updated sequentially. Generally the function consists of 7 steps:
#        1. Propose a new value for the parameter (proposed values are in lowercase, original values are in uppercase)
#        2. Update the exposure matrix
#        3. Calculate the probability of infection
#        4. Calculate the likelihood under the proposed value
#        5. Run the Metropolis Algorithm
#        6. Accept of reject the proposed value by drawing from a uniform distribution
#        7. If the proposed value is accepted, update the values (globally) the likelihood, and    
#              exposure matrices
#       For many parameters not all seven steps are needed. Step 2 is only needed for parameters which discribe the spread of the disease agent for instance. 
#		The steps are noted in the functions below.
#
#	The following functions are included:
#       -Update R		--Rate of Spread
#       -Update B		--Baseline instantaneous probablility of Infection given exposure
#       -Update RF		--Coefficient of the risk factor describing its multiplicative effect on Baseline Risk. 
#			Note: Alternative functions are include for cases with a single risk factor and multiple risk factors
#       -Update IntroT  --Time of Introduction for each epicenter; an array
#       -Update EPI		--Location (indicated by house code) of each epicenter; an array
#===============================================#


#====================================================================================#
#	Update function for R, the speed of spread of the parasite
#====================================================================================#
updateR<-function()

{

#1. Propose a new value, making sure it is >0

         r<-rnorm(1,R,sd=varR^.5)                    

   while(r<=0) r<-rnorm(1,R,sd=varR^.5)

#2. Update the exposure matrix, and find the exposure time for each individual

   p_exp<-EXP

   for (i in 1:epicenters)

   {

       p_exp[,i]<-(IntroT[i]-(DIS[,i]/r))

       p_exp[p_exp<0]<-0

   }

   maxexp<-get_maxEXP(p_exp)

#3. Calculate the probability of infection given the exposure time and covariates

   prob<-getprob(maxexp,B,RF)

#4. Calculate the likelihood under the proposed value of r, if the likelihood is zero skip the Metropolis Step

   p_lik<-lik(prob,maxexp)

   if (p_lik>0)

#5 Metropolis Algorithm

        {

      likelihood_ratio<-p_lik/LIK     

prior_ratio<-R_prior_prob(r)   /  R_prior_prob(R)

proposal_ratio<- R_proposal(R,r) /  R_proposal(r,R)

      accept_prob<-likelihood_ratio * prior_ratio * proposal_ratio

#6 Accept of reject the proposed value stochastically by drawing from a uniform distribution

       if (runif(1)<accept_prob)

#7 If the proposed value is accepted, update the values (globally) for R, the likelihood, and exposure matrices

       {

               R<<-r

               LIK<<-p_lik

               EXP<<-p_exp

               maxEXP<<-maxexp

       }

       }

 }

   
#====================================================================================#
#	Update function for B, Baseline instantaneous probablility of Infection given exposure
#====================================================================================#
updateB<-function()

     {

#1. Propose a new value for b

	b<-rnorm(1,B,sd=varB^.5)    

#2  Not needed- B does not affect exposure time

#3. Recalculate the probability of infection with the proposed value

	prob<-getprob(maxEXP,b,RF)
	
#4. Calculate the likelihood under the proposed value

	b_lik<-lik(prob,maxEXP)            #gets the proposed likelihood

#5. Run the Metropolis Algorithm if the proposed likelihood is > 0

    if(b_lik>0)

         {

			likelihood_ratio<-b_lik / LIK

			prior_ratio<-B_prior_prob(b) / B_prior_prob(B)

			proposal_ratio<-B_proposal(B,b) / B_proposal(b,B)

			accept_prob<-likelihood_ratio * proposal_ratio * prior_ratio   

#6. Accept of reject the proposed value statistically by drawing from a uniform distribution  

			if (runif(1)<=accept_prob)
#7 If the proposed value is accepted, update the values (globally)
				{
				B<<-b
				LIK<<-b_lik
				}
         }

   }


#====================================================================================#
#Update the coefficients on risk factors-version 1: Function for a single risk factor
#====================================================================================#
 if(num_covariates==1)

 {
    updateRF<-function()
     {
#1. Propose a new value for rf
		rf<-rnorm(1,RF,sd=varRF^.5)     
#2  Not needed- RF does not affect exposure time
#3. Recalculate the probability of infection with the proposed value
		prob<-getprob(maxEXP,B,rf)  
#4. Calculate the likelihood under the proposed value
		rf_lik<-lik(prob,maxEXP)             
#5. Run the Metropolis Algorithm if the proposed likelihood is > 0
		if(rf_lik>0)
			{
			likelihood_ratio<-rf_lik/LIK
			prior_ratio<-RF_prior_prob(rf)/RF_prior_prob(RF)
			proposal_ratio<-RF_proposal(RF,rf) / RF_proposal(rf,RF)
			accept_prob<-likelihood_ratio * proposal_ratio * prior_ratio   
#6. Accept of reject the proposed value statistically by drawing from a uniform distribution
			if(runif(1)<=accept_prob)
#7 If the proposed value is accepted, update the values (globally)
				{
				RF<<-rf
				LIK<<-rf_lik
				}
             }
         }
	}

   
#version 2: Function for multiple risk factor

 if(num_covariates>=1)
    {
    updateRF<-function()
     {  
       for(i in 1:num_covariates)
         {
#1. Propose a new value, rfs, for multiple covariates 
			rfs<-RF
			rfs[i]<-rnorm(1,RF[i],sd=varRF[i]^.5)   
#2  Not needed- RF does not affect exposure time
#3. Recalculate the probability of infection with the proposed value
			prob<-getprob(maxEXP,B,rfs)   
#4. Calculate the likelihood under the proposed value
			rf_lik<-lik(prob,maxEXP)             
#5. Run the Metropolis Algorithm if the proposed likelihood is > 0
			if(rf_lik>0)
				{
				likelihood_ratio<-rf_lik/LIK
				prior_ratio<-RF_prior_prob(rfs[i])/RF_prior_prob(RF[i])
				proposal_ratio<-RF_proposal(RF[i],rfs[i],varRF[i]) / RF_proposal(rfs[i],RF[i],varRF[i])
				accept_prob<-likelihood_ratio * proposal_ratio * prior_ratio   
#6. Accept of reject the proposed value statistically by drawing from a uniform distribution
				if(runif(1)<=accept_prob)
#7 If the proposed value is accepted, update the values (globally)
					{
					RF[i]<<-rfs[i]
					LIK<<-rf_lik
					}
				}
             }
		}
	}

     

 

#====================================================================================#
#	Update IntroT  --Time of Introduction for each epicenter
#		-the Introduction time is updated for each epicenter separately. 
#====================================================================================#


updateIntroT<-function(n)
    {

#1. Propose a new value, making sure it is >0
	tnew<-rnorm(1,IntroT[n],sd=varT[n]^.5)     
	while(tnew<0|tnew>40)  tnew<-rnorm(1,IntroT[n],sd=varT[n]^.5)
#2. Update the exposure matrix, and find the exposure time for each individual       
	p_exp<-EXP
	p_exp[,n]<-(tnew-(DIS[,n]/R))
	p_exp[,n][p_exp[,n]<0]<-0
	maxexp<-get_maxEXP(p_exp)
#3. Calculate the probability of infection given the exposure time and covariates             
		tnew_prob<-getprob(maxexp,B,RF)
        tnew_lik<-lik(tnew_prob,maxexp)
		if(tnew_lik>0)
			{
#5 Metropolis Algorithm
			likelihood_ratio<-tnew_lik/LIK
			prior_ratio<-T_prior_prob(tnew)/T_prior_prob(IntroT[n])
			proposal_ratio<- T_proposal(IntroT[n],tnew,varT[n]) / T_proposal(tnew,IntroT[n],varT[n])
			accept_prob<-likelihood_ratio * proposal_ratio * prior_ratio
#6. Accept of reject the proposed value statistically by drawing from a uniform distribution
            if(runif(1)<=accept_prob)
#7 If the proposed value is accepted, update the values (globally)
				{
				IntroT[n]<<-tnew
				LIK<<-tnew_lik
				EXP<<-p_exp
				maxEXP<<-maxexp
				}
			}

	}


#====================================================================================#
#	Update EPI--Location (indicated by house code) of each epicenter
#		-the location is updated for each epicenter separately. 
#====================================================================================#

updateEPI<-function(n)
    {
#1. Propose a new household as an epicenter.  The proposal distribution is a Gaussian function of the distance
#     between each candidate epicenter (each household) and the current epicenter
		proposals<-dnorm(hDIST[,epi[n]],0,sd=varEPI^.5)
		pick<-sum(proposals)*runif(1)
		newep<-which.max(cumsum(proposals)>pick)
#2a. Update the distance matrix, between each household and each epicenter.
		dis<-DIS
		dis[,n]<-dise(newep)
#2. Update the exposure matrix, given the new distance matrix, and find the exposure time for each individual           
		p_exp<-EXP
		p_exp[,n]<-(IntroT[n]-(dis[,n]/R))
		p_exp[,n][p_exp[,n]<0]<-0
		maxexp<-get_maxEXP(p_exp)
#3. Calculate the probability of infection given the exposure time and covariates                                 
		prob<-getprob(maxexp,B,RF)
#4. Calculate the likelihood under the proposed value of r, if the likelihood is zero skip the Metropolis Step           
		p_lik<-lik(prob,maxexp)
		if (p_lik>0)
#5 Metropolis-Hastings Algorithm	
			{
			likelihood_ratio<-p_lik/LIK
			prior_ratio<-EPI_prior_prob(newep)/EPI_prior_prob(epi[n])
			proposal_ratio<- EPI_proposal(epi[n],newep,varEPI) /  EPI_proposal(newep, epi[n], varEPI)
			accept_prob<-likelihood_ratio * proposal_ratio * prior_ratio   
#6. Accept of reject the proposed value statistically by drawing from a uniform distribution	
				if (runif(1)<=accept_prob)
					{
					epi[n]<<-newep	
					LIK<<-p_lik
					DIS<<-dis
					EXP<<-p_exp
					maxEXP<<-maxexp
					}
             }

    }

       

#===============================================#
#    MCMC Wrapper function
#       1. runall function to do everything given the length of the chains (reps) and the number of epicenters
#===============================================#

runall<-function(reps,epicenters)
	{         
#set the  number of epicenters in the model as a global variable   
	epicenters<<-epicenters
    reps<<-reps
#Set initial values for the parameters
    init_values(epicenters, init_B, init_RF ,init_R)
#Initialize chains
	chain_init(reps/thin,epicenters)
#Set the hyperparameters for the proposal distributions
	proposal_init(epicenters)
#Calculate the likelihood with the initial parameters, and make sure it is > 0
	lik_init(epicenters)
    while(LIK<=0)
       {
       init_values(epicenters, init_B, init_RF ,init_R)
       lik_init(epicenters)
       }
#Update the chains  
       for (i in 1:reps)
         {
			updateR()
			updateB()
			updateRF()
			for (j in 1:epicenters)
               {
				updateEPI(j)   
				updateIntroT(j)
               }
#Save values of thinned chains     
			if(i==THIN[it])
				{
				Rchain[it]<<-R
				Bchain[it]<<-B
				ifelse(num_covariates==1, RFchain[it]<<-RF, RFchain[it,]<<-RF)
				Lchain[it]<<-LIK
					for (c in 1:epicenters)
					{
					Echain[it,c]<<-epi[c]
					Tchain[it,c]<<-IntroT[c]
					}
				it<<-it+1
				}
           }
#Prepare output as a named data frame 
    if(epicenters>1)
		{
		colnames(Tchain)<-paste("T",c(1:epicenters),sep="")
		colnames(Echain)<-paste("E",c(1:epicenters),sep="")
		}
	if(num_covariates>1)
		{
		colnames(RFchain)<-paste("RF",c(1:num_covariates),sep="")
		}
	outp<-data.frame(Tchain,Echain,Rchain,Bchain,RFchain,Lchain)    
	return(outp)        
 }

  
CHAINS1<-runall(reps,epicenters)

#==============================================
# Run accompanying movie_simulated_data.r
# to visualize the MCMC fitting process.
#==============================================


