//#include <ctime>
//#include <iomanip>
//#include <iostream>
//#include "constants.h"
//#include "Mixture.h"
//#include "StifflSolver.h"
//using namespace std;
//
//int main()
//{
//	Mixture *mix = new Mixture("Substances_H2_O2_O3.txt", 
//		"Reactions_H2_O2_O3_19.txt", 
//		"MoleFractions_H2_O2_O3.txt");
//	
//	cout << "Substances" << endl;
//	for (int i = 0; i < mix->getNSpecies(); ++i) {
//		cout << setw(3) << i+1 << " " << mix->substances[i]->nameOfSubstance;
//		cout << endl;
//	}
//	cout << endl;
//
//	cout << "Reactions" << endl;
//	for (int i = 0; i < mix->nReactions; ++i) {
//		cout << setw(3) << i+1 << " " 
//			 << mix->reactions[i].nameOfReaction << endl;
//	}
//	cout << endl;
//
//	double *initialValues = new double [mix->getNSpecies() + 1];
//	for (int i = 0; i < mix->getNSpecies() + 1; ++i) {
//		initialValues[i] = 0.0;
//	}
//
//	StifflSolver kinetics(
//		mix->getNSpecies() + 1,
//		initialValues,
//		0.0,
//		0.0,
//		1.0e-13);
//
//	kinetics.Seteps(1e-2);
//	RealType *vf = new RealType[mix->getNSpecies()];
//	mix->returnMoleFractions(vf);
//
//	clock_t start;
//	clock_t finish;
//	double workingTime;
//
//	start = clock();
//	mix->setStateWithTPX(2000, 30 * ONE_ATM, vf);
//	cout << "Results" << endl;
//	cout << "oldT = " << mix->getOldTemperature() << endl;
//	cout << "oldP = " << mix->calculateOldPressure() << endl;
//	kinetics.performIntegration(*mix, 1.0);
//	mix->updateMoleFractions(vf);
//
//	cout << "T = " << mix->getTemperature() << endl;
//	cout << "P = " << mix->calculatePressure() << endl;
//	cout << "X(H2O)  = " << vf[2]  << endl;
//	cout << "X(O2) = " << vf[6] << endl;
//	cout << "X(H2) = " << vf[1] << endl;
//	finish = clock();
//
//	workingTime = (double) (finish - start) / CLOCKS_PER_SEC;
//
//	cout << "Calculations done in " << workingTime << " s." << endl << endl;
//
//	delete mix;
//	delete [] vf;
//
//	exit(0);
//}