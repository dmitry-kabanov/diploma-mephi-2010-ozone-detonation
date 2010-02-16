/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version $Id$
*
* @section DESCRIPTION
*
* Реализует функции класса Config.
*/
#include "Config.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

Config::Config(const char *filename)
{
	readFileOfConfig(filename);
}

Config::~Config()
{
}

void Config::readFileOfConfig(const char *filename)
{
	ifstream iFile(filename);
	if (!iFile) {
		cerr << "Bad filename of config file." << endl;
		exit(-1);
	}

	// Хранит в себе целую строчку из файла.
	string str;
	// Имя параметра из конфигурационного файла.
	string parName;

	while (getline(iFile, str)) {
		if (str.find("//") == 0) {
			// Считанная строка является строкой комментария.
			continue;
		}
		
		istringstream is(str);
		is >> parName;
		
		if (parName == "dx") {
			is >> dx_;
		}
		else if (parName == "dt") {
			is >> dt_;
		}
		else if (parName == "pFront") {
			is >> pFront_;
		}
		else if (parName == "pInitial") {
			is >> pInitial_;
		}
		else if (parName == "rhoFront") {
			is >> rhoFront_;
		}
		else if (parName == "rhoInitial") {
			is >> rhoInitial_;
		}
		else if (parName == "uFront") {
			is >> uFront_;
		}
		else if (parName == "uInitial") {
			is >> uInitial_;
		}
		else if (parName == "pistonVelocity") {
			is >> pistonVelocity_;
		}
		else if (parName == "gammaInsideFront") {
			is >> gammaInsideFront_;
		}
		else if (parName == "gammaAheadFront") {
			is >> gammaAheadFront_;
		}
		else if (parName == "gammaBehindFront") {
			is >> gammaBehindFront_;
		}
		else if (parName == "cellWidthCoeff") {
			is >> cellWidthCoeff_;
		}
		else if (parName == "meshSize") {
			is >> meshSize_;
		}
		else if (parName == "initialShockWaveSize") {
			is >> initialShockWaveSize_;
		}
		else if (parName == "timeSteps") {
			is >> timeSteps_;
		}
		else if (parName == "timeStepForOutput") {
			is >> timeStepForOutput_;
		}
		else if (parName == "start") {
			is >> start_;
		}
		else if (parName == "fileOfSubstances") {
			is >> fileOfSubstances_;
		}
		else if (parName == "fileOfReactions") {
			is >> fileOfReactions_;
		}
		else if (parName == "fileOfFractions") {
			is >> fileOfFractions_;
		}
		else if (parName == "fileOfPiston") {
			is >> fileOfPiston_;
		}
		else {
			cerr << "'" << parName << "' is not a known parameter." << endl;
		}
	}
}