/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version $Id$
*
* @section DESCRIPTION
*
* Реализует функции класса Output.
*/
#include "Output.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Mixture.h"
#include "Substance.h"
using namespace std;

Output::Output(const Mixture &mix, const std::string &format, std::string path)
:	names_(8 + mix.nSubstances),
	format_(format),
	path_(path),
	nColumnsOfSpecies_(mix.getNSpecies())
{
	names_[0] = "Cell No";
	names_[1] = "Right Bound (m)";
	names_[2] = "Coordinate (m)";
	names_[3] = "Pressure (Pa)";
	names_[4] = "Velocity (m/s)";
	names_[5] = "Density (kg m-3)";
	names_[6] = "Full energy (J kg-1)";
	names_[7] = "Internal Energy (J kg-1)";
	
	for (int i = 0; i < mix.nSubstances; i++) {
		names_[8+i] = "X(" + mix.substances[i]->nameOfSubstance + ")";
	}

	if (format_ == "csv") {
		delimiter_ = ";";
	}
	else {
		cerr << "Unknown format (" << format_ << "of output file." << endl;
		exit(-1);
	}
}

Output::~Output()
{
}

void Output::plotData(int timeStep, GodunovKolganMethod &gkm)
{
	int num_digits;
	char time[32];
	// TODO: добавить проверку на существование каталога.
	char filename[] = "\\data_";
	string fullname(path_);

	fullname += filename;
	num_digits = sprintf_s(time, "%d", timeStep);
	fullname += time;
	fullname += "." + format_;

	ofstream outFile(fullname.c_str());

	if (!outFile) {
		cout << "Can't open file '" << fullname << "' for writing." << endl;
	}

	outFile.setf(ios::fixed, ios::floatfield);
	outFile.precision(9);

	writeDataLabels(outFile);

	outFile << delimiter_ << gkm.getX()[0] << endl;
	for (int i = 1; i < gkm.getMeshSize(); i++) {
		outFile << gkm.getCellNumbers()[i] << delimiter_ <<
			gkm.getX()[i]                  << delimiter_ <<
			gkm.getXCenter()[i]            << delimiter_ << 
			gkm.getP()[i]                  << delimiter_ << 
			gkm.getU()[i]                  << delimiter_ << 
			gkm.getRho()[i]                << delimiter_ <<
			gkm.getFullEnergy()[i]         << delimiter_ << 
			gkm.getIntEnergy()[i]          << delimiter_;
			for (int j = 0; j < nColumnsOfSpecies_; ++j) {
				outFile << gkm.getMoleFractions()[i][j];
				outFile << (j < (nColumnsOfSpecies_ - 1) ? delimiter_ : "");
			}
			outFile << endl;
		if (gkm.getShockWaveFront()[i] == true && 
			gkm.getShockWaveFront()[i+1] == false) {
			break;
		}
	}
	outFile.flush();
}

void Output::writeDataLabels(ofstream &outFile)
{
	for (size_t i = 0; i < names_.size(); i++) {
		outFile << names_[i] << delimiter_;
	}
	outFile << endl;
}