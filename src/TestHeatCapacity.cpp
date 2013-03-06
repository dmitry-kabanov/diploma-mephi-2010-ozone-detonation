//#include <fstream>
//#include <iostream>
//#include <string>
//#include "Mixture.h"
//#include "Substance.h"
//using namespace std;
//
//int main(int argc, char *argv[])
//{
//	Mixture *mixture = new Mixture("Substances.txt", 
//		"Reactions.txt",
//		"MoleFractions.txt");
//
//	RealType *mf = new RealType[mixture->getNSpecies()];
//	for (int i = 0; i < mixture->getNSpecies(); ++i) {
//		mf[i] = 0.0;
//	}
//
//	ofstream out("thermo.csv");
//	string delimiter(";");
//	out << "T (K)" << delimiter;
//	for (int j = 0; j < mixture->getNSpecies(); ++j) {
//		out << mixture->substances[j]->nameOfSubstance << " Cp" << delimiter;
//		out << mixture->substances[j]->nameOfSubstance << " H" << delimiter;
//	}
//	out << endl;
//
//	for (int i = 298; i <= 6000; ++i) {
//		out << i << delimiter;
//		for (int j = 0; j < mixture->getNSpecies(); ++j) {
//			out << mixture->calculateSubstanceCp(j, i) << delimiter;
//			out << mixture->calculateSubstanceEnthalpy(j, i) << delimiter;
//		}
//		out << endl;
//	}
//
//	return 0;
//}