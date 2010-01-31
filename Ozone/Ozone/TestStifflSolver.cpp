//#include <ctime>
//#include <iostream>
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
//	RealType vf[3];
//	vf[2] = 100.0;
//	vf[1] = 0.0;
//	vf[0] = 0.0;
//	clock_t start;
//	clock_t finish;
//	double workingTime;
//
//	start = clock();
//	kinetics.getMixture()->setStateWithTPX(NORMAL_TEMPERATURE, ONE_ATM, vf);
//	cout << "oldT = " << kinetics.getMixture()->getOldTemperature() << endl;
//	cout << "oldP = " << kinetics.getMixture()->calculateOldPressure() << endl;
//	kinetics.performIntegration(950.0);
//	kinetics.updateMoleFractions(vf);
//
//	cout << "T = " << kinetics.getMixture()->getTemperature() << endl;
//	cout << "P = " << kinetics.getMixture()->calculatePressure() << endl;
//	cout << "X(O)  = " << vf[0]  << endl;
//	cout << "X(O2) = " << vf[1] << endl;
//	cout << "X(O3) = " << vf[2] << endl;
//	finish = clock();
//
//	workingTime = (double) (finish - start) / CLOCKS_PER_SEC;
//
//	cout << "Calculations done in " << workingTime << " s." << endl << endl;
//
//	exit(0);
//}