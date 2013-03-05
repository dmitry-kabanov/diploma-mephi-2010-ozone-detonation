/**
 * @file
 * @author Кабанов Дмитрий <kabanovdmitry@gmail.com
 * @version $Id$
 *
 * @section DESCRIPTION
 *
 * Главный файл проекта Ozone.
 */
#include "main.h"
#include <ctime>
#include <iostream>
#include <string>
#include "GodunovKolganMethod.h"
using namespace std;

//int main(int argc, char *argv[])
//{
//    std::string configFile;
//    if (argc > 1) {
//        configFile = argv[1];
//    }
//    else {
//        configFile = "Config.txt";
//    }
//
//	Config config(configFile.c_str());
//
//	GodunovKolganMethod gkm(config);
//
//	clock_t start  = clock();
//	gkm.run();
//	clock_t finish = clock();
//
//	cout << "Calculations done in " <<
//		(double) (finish - start) / CLOCKS_PER_SEC << endl;
//
//	return 0;
//}