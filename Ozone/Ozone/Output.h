/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version $Id$
*
* @section DESCRIPTION
*
* Объявляет класс Output.
*/
#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <vector>
#include <string>
#include "GodunovKolganMethod.h"
#include "Mixture.h"

/**
 * Вывод данных в файл.
 */
class Output
{
public:
	/**
	 * Конструктор класса.
	 *
	 * @param mix		газовая смесь
	 * @param format	формат выходного файла
	 * @param path		путь к каталогу выходных файлов
	 */
	Output(const Mixture &mix, const std::string &format, std::string path);
	/**
	 * Деструктор класса.
	 */
	~Output();
	/**
	 * Выводит данные в файл.
	 *
	 * @param timeStep	временной шаг 
	 * @param gkm		объект, данные которого выводятся в файл
	 */
	void plotData(int timeStep, GodunovKolganMethod &gkm);
	/**
	 * Выводит заголовки столбцов данных в файл.
	 *
	 * @param outFile выходной файл.
	 */
	void writeDataLabels(std::ofstream &outFile);

private:
	std::vector<std::string> names_;
	std::string format_;
	std::string delimiter_;
	std::string path_;
};

#endif // OUTPUT_H
