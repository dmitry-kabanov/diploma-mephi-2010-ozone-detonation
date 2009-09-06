#include "ReactionDataVar.h"
#include <cmath>

#define Rgas_cal 0.001985846 //8.31434/4.1868/1E3

 ReactionDataVar::ReactionDataVar()
{

    n=0;
	  Tlow =NULL;
	  Tup = NULL;
	  N = NULL;
	  LogA = NULL;
	  Ea = NULL;
	 
	  LogA_low = NULL;
	  N_low = NULL;
	  Ea_low = NULL;
	  LogA_Hi = NULL;
		N_Hi = NULL;
		Ea_Hi = NULL; 
		Tr1 = NULL;
		Tr2 = NULL;
		Tr3 = NULL; 
		Tr4 = NULL;
		Sri1 = NULL;
		Sri2 = NULL;
		Sri3 = NULL; 
		Sri4 = NULL; 
		Sri5 =NULL;
		bM_Type =false;
		pEff=NULL;
		nNumberOfColliderInReagentList=0;
}

ReactionDataVar::~ReactionDataVar()
{
	if(n>0)
	{
		delete [] Tlow;
		Tlow =NULL;
		delete [] Tup;
		Tup=NULL;
		delete [] N;
		N=NULL;
		delete [] LogA;
		LogA=NULL;
		delete [] Ea;
		Ea=NULL;
	
		delete [] LogA_low;
		delete [] N_low; 
		delete [] Ea_low;
		delete [] LogA_Hi;
		delete [] N_Hi;
		delete [] Ea_Hi;
		delete [] Tr1;
		delete [] Tr2;
		delete [] Tr3; 
		delete [] Tr4;
		delete [] Sri1;
		delete [] Sri2;
		delete [] Sri3;
		delete [] Sri4;
		delete [] Sri5;

		LogA_low= NULL;
		N_low= NULL; 
		Ea_low= NULL;
		LogA_Hi= NULL;
		N_Hi= NULL;
		Ea_Hi= NULL;
		Tr1= NULL;
		Tr2= NULL;
		Tr3= NULL; 
		Tr4= NULL;
		Sri1= NULL;
		Sri2= NULL;
		Sri3= NULL;
		Sri4= NULL;
		Sri5= NULL;

		delete [] pEff;
		pEff=NULL;
	}
}

int ReactionDataVar::AllocateMemoryForTemperatureRange()
{
		Tlow = new double[n];
		Tup = new double[n];
		Ea= new double[n];
		N  = new double[n];
		LogA = new double[n];
		return 0;
}
int ReactionDataVar::AllocateMemoryForCollission(int nSp)
{
		pEff = new double[nSp];
		for(int i=0;i<nSp;i++)
			pEff[i]=1.0;
		return 0;
}
int ReactionDataVar::AllocateMemoryForTemperatureRange_P()
{
		Tlow = new double[n];
		Tup = new double[n];	
		LogA_low=new double[n];
		N_low=new double[n];
		Ea_low=new double[n];
		LogA_Hi=new double[n];
		N_Hi=new double[n];
		Ea_Hi=new double[n];
		Tr1=new double[n];
		Tr2=new double[n];
		Tr3=new double[n];
		Tr4=new double[n];
		Sri1=new double[n];
		Sri2=new double[n];
		Sri3=new double[n];
		Sri4=new double[n];
		Sri5= new double[n];
		
		return 0;
}
double ReactionDataVar::GetReactionConstant(double T, double sumConc)
{
	double vel;

	double sum_gibbs=0;

	int j;
	if(n<1) return 0.0;
	if(1==n) j=0;
	else if(T>Tup[n-1]) j=n-1;
	else
	{
		for(j=n-1; j>=0;j--)
		{
			if(	T>=Tlow[j]) break;
		}
		if(-1==j)
			j=0;		
	}
	
	if(0==nTypeOFReaction)
	{
			
		vel=exp(logl(T) * N[j] -Ea[j]/Rgas_cal/T + 2.30258*LogA[j]);
	
		if(bM_Type)
		
		vel*=sumConc;
	}
	else if(11<=nTypeOFReaction  && 12>=nTypeOFReaction)
	{	

		
		
		double rate_inf;
		double rate_low;
		
		double P;

		rate_low= exp(logl(T)*N_low[j]-Ea_low[j]/(T*Rgas_cal) +2.30258*LogA_low[j]);

		rate_inf= exp(logl(T)*N_Hi[j]-Ea_Hi[j]/(T*Rgas_cal) +2.30258*LogA_Hi[j]);


		P=rate_low*sumConc/rate_inf;

		vel=(rate_low)/(1+P);

		if(bM_Type )
			vel*=sumConc;
	
		if(12==nTypeOFReaction)
		{	
			double Fcent=(1-Tr1[j])*exp(-T/Tr2[j])
								+Tr1[j]*exp(-T/Tr3[j])
								 +exp(-(Tr4[j])/T);
			double C=-0.4-0.67*log(Fcent);
			double N=0.75-1.27*log(Fcent);
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			P=1.35344029177355E-02 ;
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!
			double F=exp((1/(1+(log(P)+C)*(log(P)+C)
			/(N-0.14*(log(P)+C))/(N-0.14*(log(P)+C)))*log(Fcent)));

			vel*=F;
		}
			if(13==nTypeOFReaction)
		{	
			double x = 1/(1+log(P*P)); // check !!!!!
			double F = pow(
				Sri1[j]*exp(-Sri2[j]/T) +exp(-T/Sri3[j]),-x)
				*Sri4[j]*pow(T,Sri5[j]);

			vel*=F;
		}
	}
		return vel;
}


