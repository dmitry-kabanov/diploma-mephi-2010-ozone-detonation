//#include <ctime>
//#include <iostream>
//#include "constants.h"
//#include "StifflSolver.h"
//using namespace std;
//
//int main()
//{
//	Mixture *mix = new Mixture("Substances_O3.txt", "Reactions_O3.txt",
//		"MoleFractions_O3.txt");
//	for (int i = 0; i < mix->getNSpecies(); ++i) {
//		cout << mix->reactions[i].nameOfReaction << endl;
//	}
//
//	double *vf1 = new double[mix->getNSpecies() + 1];
//	for (int i = 0; i < mix->getNSpecies() + 1; ++i) {
//		vf1[i] = 0.0;
//	}
//	StifflSolver kinetics(
//		mix->getNSpecies() + 1,
//		vf1,
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
//	//kinetics.getMixture()->setStateWithTPX(1909, 6543160, vf);
//	//mix->setStateWithURhoX(4445030, 19.7335, vf);
//	mix->setStateWithURhoX(4482555, 19.7335, vf);
//	double temperature = 1914;
//	cout << "Mixture Enthalpy at T = " << temperature <<
//		    " equals to " << mix->calculateMixtureEnthalpy(temperature) << endl;
//	cout << "oldT = " << mix->getOldTemperature() << endl;
//	cout << "oldP = " << mix->calculateOldPressure() << endl;
//	kinetics.performIntegration(*mix, 1.0);
//	mix->updateMoleFractions(vf);
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