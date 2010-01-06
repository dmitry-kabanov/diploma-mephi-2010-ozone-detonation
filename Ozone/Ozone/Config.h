/**
* @file
*
* @author  ������� ������� <kabanovdmitry@gmail.com>
* @version $Id$
*
* @section DESCRIPTION
*
* ��������� ����� Config.
*/
#ifndef CONFIG_H
#define CONFIG_H

#include "RealType.h"
#include <string>

class Config
{
public:
	/**
	 * �����������.
	 *
	 * @param filename ��� ����������������� �����.
	 */
	Config(const char *filename);
	/**
	 * ����������.
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
	bool getResume() const { return resume_; }
	void setResume(bool val) { resume_ = val; }
	int getStart() const { return start_; }
	void setStart(int val) { start_ = val; }
	std::string getFileOfSubstances() const { return fileOfSubstances_; }
	void setFileOfSubstances(std::string val) { fileOfSubstances_ = val; }
	std::string getFileOfReactions() const { return fileOfReactions_; }
	void setFileOfReactions(std::string val) { fileOfReactions_ = val; }
	std::string getFileOfFractions() const { return fileOfFractions_; }
	void setFileOfFractions(std::string val) { fileOfFractions_ = val; }

private:
	/**
	 * ��������� ������������ �� ��������� �����.
	 *
	 * @param filename ��� ����������������� �����.
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

	// ������ ���������� �����, �����.
	int meshSize_;
	// ����� �����, ���������� ������� ����� � ��������� ������ �������.
	int initialShockWaveSize_;
	// ������ ���������� ����� �� ������� � �������.
	int timeSteps_;
	// ���������� ��������� �����, ����� ������� ���������� ������������ � ����.
	int timeStepForOutput_;
	// ����, ���������� ��� ������ ������ ������. ������ � ���������� start, 
	// ������� ������ ��� �������, �� ������� ���������� ������.
	// ���� resume ����� ����, �� start ������ ���� ����� 1.
	// ���� resume ����� �������, ��, ��������, start, ������ 5000, ��������, ���
	// ��������� �������� ������ �� �����, ����������� 5000-� ��������� ���,
	// � ������ ������ � 5001-�� ����.
	bool resume_;
	int start_;
	// ��� �����, ����������� ������ � ���������.
	std::string fileOfSubstances_;
	// ��� �����, ����������� ������ � ��������.
	std::string fileOfReactions_;
	// ��� �����, ����������� ������� ���� ����������� ��������� �����.
	std::string fileOfFractions_;
};

#endif // CONFIG_H
