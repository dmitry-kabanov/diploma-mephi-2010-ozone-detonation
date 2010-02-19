/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version $Id$
*
* @section DESCRIPTION
*
* Объявляет класс Config.
*/
#ifndef CONFIG_H
#define CONFIG_H

#include "RealType.h"
#include <string>

class Config
{
public:
	/**
	 * Конструктор.
	 *
	 * @param filename имя конфигурационного файла.
	 */
	Config(const char *filename);
	/**
	 * Деструктор.
	 */
	~Config();

	RealType getDx() const { return dx_; }
	void setDx(RealType val) { dx_ = val; }
	RealType getDt() const { return dt_; }
	void setDt(RealType val) { dt_ = val; }
	RealType getPFront() const { return pFront_; }
	void setPFront(RealType val) { pFront_ = val; }
	RealType getPInitial() const { return pInitial_; }
	void setPInitial(RealType val) { pInitial_ = val; }
	RealType getRhoFront() const { return rhoFront_; }
	void setRhoFront(RealType val) { rhoFront_ = val; }
	RealType getRhoInitial() const { return rhoInitial_; }
	void setRhoInitial(RealType val) { rhoInitial_ = val; }
	RealType getUFront() const { return uFront_; }
	void setUFront(RealType val) { uFront_ = val; }
	RealType getUInitial() const { return uInitial_; }
	void setUInitial(RealType val) { uInitial_ = val; }
	RealType getPistonVelocity() const { return pistonVelocity_; }
	void setPistonVelocity(RealType val) { pistonVelocity_ = val; }
	RealType getGammaInsideFront() const { return gammaInsideFront_; }
	void setGammaInsideFront(RealType val) { gammaInsideFront_ = val; }
	RealType getGammaAheadFront() const { return gammaAheadFront_; }
	void setGammaAheadFront(RealType val) { gammaAheadFront_ = val; }
	RealType getGammaBehindFront() const { return gammaBehindFront_; }
	void setGammaBehindFront(RealType val) { gammaBehindFront_ = val; }
	RealType getCellWidthCoeff() const { return cellWidthCoeff_; }
	void setCellWidthCoeff(RealType val) { cellWidthCoeff_ = val; }
	int getMeshSize() const { return meshSize_; }
	void setMeshSize(int val) { meshSize_ = val; }
	int getInitialShockWaveSize() const { return initialShockWaveSize_; }
	void setInitialShockWaveSize(int val) { initialShockWaveSize_ = val; }
	int getTimeSteps() const { return timeSteps_; }
	void setTimeSteps(int val) { timeSteps_ = val; }
	int getTimeStepForOutput() const { return timeStepForOutput_; }
	void setTimeStepForOutput(int val) { timeStepForOutput_ = val; }
	int getStart() const { return start_; }
	void setStart(int val) { start_ = val; }
	std::string getFileOfSubstances() const { return fileOfSubstances_; }
	void setFileOfSubstances(std::string val) { fileOfSubstances_ = val; }
	std::string getFileOfReactions() const { return fileOfReactions_; }
	void setFileOfReactions(std::string val) { fileOfReactions_ = val; }
	std::string getFileOfFractions() const { return fileOfFractions_; }
	void setFileOfFractions(std::string val) { fileOfFractions_ = val; }
	std::string getFileOfPiston() const { return fileOfPiston_; }
	void setFileOfPiston(std::string val) { fileOfPiston_ = val; }
	RealType getMinConcentration() const { return minConcentration_; }
	void setMinConcentration(RealType val) { minConcentration_ = val; }
	int getWhatSpecies() const { return whatSpecies_; }
	void setWhatSpecies(int val) { whatSpecies_ = val; }
	int getNCellsBehindLeadShock() const { return nCellsBehindLeadShock_; }
	void setNCellsBehindLeadShock(int val) { nCellsBehindLeadShock_ = val; }
	int getNCellsInReductionZone() const { return nCellsInReductionZone_; }
	void setNCellsInReductionZone(int val) { nCellsInReductionZone_ = val; }

private:
	/**
	 * Считывает конфигурацию из заданного файла.
	 *
	 * @param filename имя конфигурационного файла.
	 */
	void readFileOfConfig(const char *filename);

	RealType dx_;
	RealType dt_;
	RealType pFront_;
	RealType pInitial_;
	RealType rhoFront_;
	RealType rhoInitial_;
	RealType uFront_;
	RealType uInitial_;

	RealType pistonVelocity_;

	RealType gammaInsideFront_;
	RealType gammaAheadFront_;
	RealType gammaBehindFront_;

	RealType cellWidthCoeff_;

	// Размер разностной сетки, ячеек.
	int meshSize_;
	// Число ячеек, содержащих ударную волну в начальный момент времени.
	int initialShockWaveSize_;
	// Полное количество шагов по времени в расчете.
	int timeSteps_;
	// Количество временных шагов, через которые результаты записываются в файл.
	int timeStepForOutput_;
	// Временной шаг, с которого необходимо начать интегрирование.
	// Если start_ равно 0, то расчет начинается с начала. 
	// Если start_, например, равно 5000, то программа загрузит данные 
	// из файла с именем "data_5000.txt", содержащего 5000-й временной шаг,
	// и начнет расчет с 5001-го шага.
	int start_;
	// Имя файла, содержащего данные о веществах.
	std::string fileOfSubstances_;
	// Имя файла, содержащего данные о реакциях.
	std::string fileOfReactions_;
	// Имя файла, содержащего мольные доли компонентов начальной смеси.
	std::string fileOfFractions_;
	// Имя файла, содержащего данные о поршне.
	std::string fileOfPiston_;

	// Концентрация вещества в зоне объединения, при которой объединение
	// можно проводить.
	RealType minConcentration_;
	// Номер вещества, по которому смотрится концентрация.
	int whatSpecies_;
	// Количество ячеек за фронтом лидирующей ударной волны, в которых
	// объединение ячеек производить нельзя.
	int nCellsBehindLeadShock_;
	// Минимальное количество ячеек в зоне их объединения,
	// при котором объединение можно проводить.
	int nCellsInReductionZone_;
};

#endif // CONFIG_H
