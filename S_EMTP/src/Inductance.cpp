#include "Inductance.h"
#include <iostream>
using namespace std;

Inductance::Inductance(int id,int fromNode,int toNode,double value)
{
	type=2;
	this->id=id;
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	inductanceValue=value;
	nortonEquivalentResistance = 2*inductanceValue/deltaT;
	need_NEC=1;
	isSeriesOrNot = isSeries();
}

void Inductance::calculateNortonEquivalentCurrent(double time)
{//计算支路的诺顿等效电路中的电流项
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
	nortonEquivalentCurrent = branchCurrent + branchVoltage/nortonEquivalentResistance;
}

void Inductance::calculateNortonEquivalentResistance(double time)
{//计算支路的诺顿等效电阻
	nortonEquivalentResistance = 2*inductanceValue/deltaT;
}