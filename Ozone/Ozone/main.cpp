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
#include "GodunovKolganMethod.h"

int main(/*int argc, char *argv[]*/)
{
	Config config("Config.txt");

	GodunovKolganMethod gkm(config);
	gkm.run();

	return 0;
}