#include "Mixture.h"
#include <cmath>
#include <fstream>

Mixture::Mixture()
{
	nameOfFileOfErrors = NULL;
}

Mixture::~Mixture()
{
	delete [] Gibbs_Ef;
	delete [] MU;
	delete [] reactionsDatVar;
	delete [] H0;
	delete [] thermDatVar;
	delete [] namesOfSubstances;
	delete [] Atom;
}


int Mixture::readFileOfSubstances(const char *fileName)
{	
	char b0;
	int numberOfSubstances;

    std::ifstream File_Of_Substances(fileName, 
                                     std::ios::in);
	if(0 == File_Of_Substances.is_open())
	{
        std::string str2;
		str2 = "Can not open file ";
        str2 += fileName;
        strcpy_s(ErrorBuffer, str2.c_str());
		return 11;
	}

	File_Of_Substances >> numberOfSubstances >> b0;

	nSpecies = numberOfSubstances;

	thermDatVar          = new ThermDataVar[nSpecies];
	Gibbs_Ef             = new double[nSpecies];
	MU                   = new double[nSpecies];
	H0                   = new double[nSpecies];
    namesOfSubstances    = new std::string[nSpecies];
	Atom                 = new Elements[nSpecies];
	
	int nLastElement = 0, nCharNumber = -1, nCharElementsNumber = -1, nKoeff = 0;	
	
    for(int i =0;i<numberOfSubstances;i++)
	{
		nLastElement=0, nCharNumber=-1, nCharElementsNumber=-2;
		File_Of_Substances>>b0;	
	
		while(b0!=' ')
		{
			nCharNumber++;
			namesOfSubstances[i]+=b0;
					if(nCharNumber==nCharElementsNumber+1) nKoeff=1;
					else nKoeff=0;
							 
					if(b0=='0') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10;
					else if(b0=='1') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10+1;
					else if(b0=='2') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10+2;
					else if(b0=='3') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10+3;
					else if(b0=='4') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10+4;
					else if(b0=='5') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10+5;
					else if(b0=='6') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10+6;
					else if(b0=='7') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10+7;
					else if(b0=='8') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10+7;
					else if(b0=='9') Atom[i].nBeta[nLastElement] = (Atom[i].nBeta[nLastElement]-nKoeff)*10+9;
					else
					{
						
						for (int j=0;j<Atom[i].N;j++)
						{
							if(b0==Atom[i].ElementsList[j]) 
							{
 								nLastElement=j;
								nCharElementsNumber=nCharNumber;
								Atom[i].nBeta[nLastElement]=1;
							}
						}
					}
			 		
			File_Of_Substances.read(&b0,1);;
		}
		
			File_Of_Substances>>b0;
			File_Of_Substances>>b0;
	

		File_Of_Substances>>H0[i]>>MU[i];

		File_Of_Substances>>thermDatVar[i].n;

		thermDatVar[i].AllocateMemoryForTemperatureRange();

		for(int j=0;j<thermDatVar[i].n;j++)
		{
			File_Of_Substances>>thermDatVar[i].Tlow[j];
			File_Of_Substances>>thermDatVar[i].Tup[j];
			for(int j1=0;j1<7;j1++)
				File_Of_Substances>>thermDatVar[i].a[j][j1];
		}

			File_Of_Substances>>b0; //0
			File_Of_Substances>>b0; //0
			for(int j1=0;j1<4;j1++)
				File_Of_Substances>>b0; //0 //USER

	}
	File_Of_Substances.close();
	return 0;
}

int Mixture::readFileOfReactions(const char *fileName)
{

	int  NumberOfReactions;
	char b0;

    std::ifstream File_Of_Reactions(fileName, std::ios::in);
	
	if (0 == File_Of_Reactions.is_open())
	{
        std::string str2;
		str2 = "Can not open file ";
        str2 += fileName;
        strcpy_s(ErrorBuffer, str2.c_str());
		return 11;
	}

	File_Of_Reactions>>NumberOfReactions>>b0;

	nReactions = NumberOfReactions;

	reactionsDatVar = new ReactionDataVar[nReactions+1];
		
	for(int i=0;i<NumberOfReactions;i++)
	{
		File_Of_Reactions>>reactionsDatVar[i].nTypeOFReaction;
		if(0==reactionsDatVar[i].nTypeOFReaction)
		{
			File_Of_Reactions>>reactionsDatVar[i].nNumberOfColliderInReagentList; // take another 0

			File_Of_Reactions>>reactionsDatVar[i].n;

			File_Of_Reactions>>b0; //Take1 1 it is just need
			if(reactionsDatVar[i].n>0)
				reactionsDatVar[i].AllocateMemoryForTemperatureRange();

			for(int j=0;j<reactionsDatVar[i].n;j++)
			{
				File_Of_Reactions>>reactionsDatVar[i].Tlow[j];
				File_Of_Reactions>>reactionsDatVar[i].Tup[j];
				File_Of_Reactions>>reactionsDatVar[i].LogA[j];
				File_Of_Reactions>>reactionsDatVar[i].N[j];
				File_Of_Reactions>>reactionsDatVar[i].Ea[j];

				
			}

			File_Of_Reactions>>b0;
			while(b0!=10)
			{
				
				reactionsDatVar[i].NameOfReaction+=b0;
				File_Of_Reactions.read(&b0,1);
			}
		}
	
		////pressure dependence reaction
		//if(1==reactionsDatVar[i].nTypeOFReaction)
		//{
		//	File_Of_Reactions>>reactionsDatVar[i].nNumberOfColliderInReagentList;

		//	File_Of_Reactions>>reactionsDatVar[i].n;

		//	File_Of_Reactions>>b0; //Take 1 it is just need
		//	if(reactionsDatVar[i].n>0)
		//		reactionsDatVar[i].AllocateMemoryForTemperatureRange_P();


		//	for(int j=0;j<reactionsDatVar[i].n;j++)
		//	{
		//		File_Of_Reactions>>
		//		reactionsDatVar[i].Tlow[j]>>
		//		reactionsDatVar[i].Tup[j]>>
		//		reactionsDatVar[i].LogA_low[j]>>
		//		reactionsDatVar[i].N_low[j]>>
		//		reactionsDatVar[i].Ea_low[j]>>
		//		reactionsDatVar[i].LogA_Hi[j]>>
		//		reactionsDatVar[i].N_Hi[j]>>
		//		reactionsDatVar[i].Ea_Hi[j]>>
		//		reactionsDatVar[i].Tr1[j]>>
		//		reactionsDatVar[i].Tr2[j]>>
		//		reactionsDatVar[i].Tr3[j]>>
		//		reactionsDatVar[i].Tr4[j]>>
		//		reactionsDatVar[i].Sri1[j]>>
		//		reactionsDatVar[i].Sri2[j]>>
		//		reactionsDatVar[i].Sri3[j]>>
		//		reactionsDatVar[i].Sri4[j]>>
		//		reactionsDatVar[i].Sri5[j]; 
		//	}

		//	File_Of_Reactions>>b0;
		//	while(b0!=10)
		//	{
		//		
		//		reactionsDatVar[i].NameOfReaction+=b0;
		//		File_Of_Reactions.read(&b0,1);
		//	}
		//	
		//	
		//	if(reactionsDatVar[i].Tr1[0]>= 1.0E-50) reactionsDatVar[i].nTypeOFReaction=12;
		//	else if(reactionsDatVar[i].Sri1[0]>= 1.0E-50) reactionsDatVar[i].nTypeOFReaction=13;
		//	else reactionsDatVar[i].nTypeOFReaction=11;
		//}
	}
			
	File_Of_Reactions.close();
	//int nCode=FillProductsReagentsArrays();
	//if(nCode) return nCode;
	//return CheckForDoubleReactions();
}
