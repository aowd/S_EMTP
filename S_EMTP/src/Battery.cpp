#include "Battery.h"
#include <cmath>
#include <iostream>
using namespace std;

Battery::Battery(int id,int fromNode,int toNode,double initSOC)
{
	type=28;
	this->id=id;
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	this->initSOC = initSOC;
	need_NEC=1;
	isSeriesOrNot = isSeries();
	//蓄电池模型取值
	Qfull = 9600;
	Efull = 931.2;
	Eexp = 864.3;
	Qexp = 471.65;
	Enom = 800;
	Qnom = 8681.7;
	Inom = 4173.9;
	innerResistance = 0.000833;

	Qit = Qfull*(1-initSOC);
	this->innerResistance = innerResistance;
	nortonEquivalentResistance = innerResistance;
}

void Battery::calculateNortonEquivalentCurrent(double time)
{
	calculateParameters();
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
	nortonEquivalentCurrent = (E0-K*Qfull/(Qfull-Qit)+A*exp(-B*Qit))/nortonEquivalentResistance;
}

void Battery::calculateNortonEquivalentResistance(double time)
{
	nortonEquivalentResistance = innerResistance;
}

void Battery::calculateBranchCurrent()
{
	branchCurrent_1 = branchCurrent;
	branchCurrent = nortonEquivalentCurrent + branchVoltage/nortonEquivalentResistance;
	Qit+=branchCurrent*deltaT/3600;//单位为Ah
}

void Battery::calculateParameters()
{
	A = Efull-Eexp;
	B = 3/Qexp;
	K = (Efull-Enom+A*(exp(-B*Qnom)-1))*(Qfull-Qnom)/Qnom;
	E0 = Efull+K+nortonEquivalentResistance*branchCurrent-A;
}

void Battery::storeInternalVariables()
{
	branchVoltage_bak = branchVoltage;
	branchCurrent_bak = branchCurrent;
	nortonEquivalentCurrent_bak = nortonEquivalentCurrent;
	Qit_bak = Qit;
}

void Battery::restoreInternalVariables()
{
	branchVoltage = branchVoltage_bak;
	branchCurrent = branchCurrent_bak;
	nortonEquivalentCurrent = nortonEquivalentCurrent_bak;
	Qit = Qit_bak;
}