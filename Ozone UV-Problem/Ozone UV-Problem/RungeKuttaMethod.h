/**
* @file
*
* @author  ������� ������� <kabanovdmitry@gmail.com>
* @version %I%
*
* @section DESCRIPTION
*
* ����� �������������� ������� ��� ���������� ��������.
*/
#ifndef RUNGE_KUTTA_METHOD_H
#define RUNGE_KUTTA_METHOD_H

#include <iostream>
#include <fstream>
#include "RealType.h"
#include "Mixture.h"
#include "Reaction.h"
#include "Substance.h"

class RungeKuttaMethod
{
public:
    /**
     * ����������� ������.
     *
     * @param initialPressure       ��������� ��������, ��
     * @param initialTemperature    ��������� �����������, �
     * @param initialTimeStep       ��������� ��������� ���, �
     * @param aFullTime             ���������� ����� ��������������
     * @param aTimeStepForOutput    ����� ������� ��������� ����� ��������
     * @param fileOfSubstances      ��� ����� � ����������
     * @param fileOfReactions       ��� ����� � ���������
     * @param fileOfVolumeFractions ��� ����� � ��������� ������
     */
    RungeKuttaMethod(RealType    initialPressure,
                     RealType    initialTemperature,
                     RealType    initialTimeStep,
                     int         aFullTime,
                     int         aTimeStepForOutput,
                     const char *fileOfSubstances,
                     const char *fileOfReactions,
                     const char *fileOfVolumeFractions);
    /**
     * ���������� ������.
     */
    ~RungeKuttaMethod();
private:
    /**
     * ���������� �������������� ������� ���.
     */
    void performIntegration();

    /**
     * ��������� �������� �������� ������� �� O3 ��� �������� �����������
     * � �������.
     *
     * @param t         �����������, �
     * @param concOfO   ������������ O, ������� / ��**3
     * @param concOfO3  ������������ O3, ������� / ��**3
     * @param concOfO2  ������������ O2, ������� / ��**3
     * @return          �������� �������� ������� ������������ O3,
     * ������� / (��**3 * �)
     */
    RealType rightSideForO3(RealType t, 
                            RealType concOfO, 
                            RealType concOfO3,
                            RealType concOfO2);
    /**
    * ��������� �������� �������� ������� �� O ��� �������� �����������
    * � �������.
    *
    * @param t         �����������, �
    * @param concOfO   ������������ O, ������� / ��**3
    * @param concOfO3  ������������ O3, ������� / ��**3
    * @param concOfO2  ������������ O2, ������� / ��**3
    * @return          �������� �������� ������� ������������ O,
    * ������� / (��**3 * �)
    */
    RealType rightSideForO(RealType t, 
                           RealType concOfO, 
                           RealType concOfO3, 
                           RealType concOfO2);
    /**
    * ��������� �������� �������� ������� �� O2 ��� �������� �����������
    * � �������.
    *
    * @param t         �����������, �
    * @param concOfO   ������������ O, ������� / ��**3
    * @param concOfO3  ������������ O3, ������� / ��**3
    * @param concOfO2  ������������ O2, ������� / ��**3
    * @return          �������� �������� ������� ������������ O2,
    * ������� / (��**3 * �)
    */
    RealType rightSideForO2(RealType t, 
                           RealType concOfO, 
                           RealType concOfO3, 
                           RealType concOfO2);
    /**
     * ��������� ��������� �������� ������ �������.
     *
     * @param i  ���������� ����� ������� � ������� �������
     * @return   �������� ��������� �������� ������ �������.
     */
    RealType calculateRateForForwardReaction(int i);
    /**
    * ��������� ��������� �������� �������� �������.
    *
    * @param i  ���������� ����� ������� � ������� �������
    * @param kf ��������� �������� ������ �������
    * @return   �������� ��������� �������� �������� �������.
    */
    RealType calculateRateForBackReaction(int i, RealType kf);
    /**
     * ��������� ��� ��������������, �
     */
    RealType h;
    /**
     * ����� �����.
     */
    Mixture *mixture;
    /**
     * ���������� ��������� ����� ��������������.
     */
    int fullTime;
    /**
     * ���������� ��������� �����, ����� ������� ���������� ����� � ����.
     */
    int timeStepForOutput;
    /**
     * �����, �
     */
    RealType time;
    /**
     * �������� ����.
     */
    std::ofstream outputFile;
    /**
     * ������� ���������� ������� � ����.
     */
    void printToFile();
    /**
     * ������� ��������� � ����.
     */
    void printHeadingToFile();
};

#endif