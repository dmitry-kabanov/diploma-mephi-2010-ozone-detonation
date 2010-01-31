/**
 * @file
 * Это тест на то, за какое время разлагается озон.
 */
//#include <ctime>
//#include <iostream>
//#include <fstream>
//#include "constants.h"
//#include "StifflSolver.h"
//using namespace std;
//
//int main()
//{
//	int nSpecies = 3;
//	double vf1[4];
//	vf1[0] = 0.0;
//	vf1[1] = 0.0;
//	vf1[2] = 0.0;
//	vf1[3] = 0.0;
//	StifflSolver kinetics(
//		nSpecies + 1,
//		vf1,
//		0.0,
//		0.0,
//		1.0e-13,
//		"Substances.txt",
//		"Reactions.txt",
//		"MoleFractions.txt");
//
//	kinetics.Seteps(1e-2);
//	ofstream outFile("outputOzoneDecomposition.txt");
//	outFile << "t T P X(O) X(O2) X(O3)" << endl;
//	RealType vf[3];
//
//	for (int i = 1; i < 1000; i += 10) {
//		vf[2] = 100.0;
//		vf[1] = 0.0;
//		vf[0] = 0.0;
//		kinetics.getMixture()->setStateWithTPX(NORMAL_TEMPERATURE, ONE_ATM, vf);
//		kinetics.performIntegration((double) i);
//		kinetics.updateMoleFractions(vf);
//
//		outFile << i << " ";
//		outFile << kinetics.getMixture()->getTemperature() << " ";
//		outFile << kinetics.getMixture()->calculatePressure() << " ";
//		outFile << vf[0]  << " ";
//		outFile << vf[1] << " ";
//		outFile << vf[2] << endl;
//	}
//
//	exit(0);
//}