#include "Capacitance.h"
#include <iostream>
using namespace std;

Capacitance::Capacitance(int id,int fromNode,int toNode,double value)
{
	type=3;
	this->id=id;
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	capacitanceValue =value;
	nortonEquivalentResistance = deltaT/(2*capacitanceValue);
	need_NEC=1;
	isSeriesOrNot = isSeries();
}

void Capacitance::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
	nortonEquivalentCurrent = -branchCurrent - branchVoltage/nortonEquivalentResistance;
}

void Capacitance::calculateNortonEquivalentResistance(double time)
{
	nortonEquivalentResistance = deltaT/(2*capacitanceValue);
}
