#ifndef REACTION_DATA_VAR_H
#define REACTION_DATA_VAR_H

#include <string>

#define MAX_REAG 4	   //max. num. of reagents
#define MAX_PROD 5	   //max. num. of products
class ReactionDataVar  
{
public:
	int nTypeOFReaction;
	int nDirection;
	int nNumberOfColliderInReagentList;
	double *N;
	double *Ea;
	double *LogA;
	double *Tup;
	double *Tlow;
	double *pEff;


	double *LogA_low, *N_low, *Ea_low, *LogA_Hi,
		*N_Hi, *Ea_Hi, *Tr1, *Tr2, *Tr3, *Tr4, *Sri1, *Sri2,
		*Sri3, *Sri4, *Sri5;

    std::string NameOfReaction;
	int AllocateMemoryForTemperatureRange();
	int AllocateMemoryForTemperatureRange_P();
	int AllocateMemoryForCollission(int nSp);
	int n;
	ReactionDataVar();
	virtual ~ReactionDataVar();
	double GetReactionConstant(double T, double sumConc);

	int reag[MAX_REAG]; //reagents(-128 - none or third body)
	int prod[MAX_PROD];
	bool bM_Type;

};

#endif