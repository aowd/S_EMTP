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
{//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
	nortonEquivalentCurrent = branchCurrent + branchVoltage/nortonEquivalentResistance;
}

void Inductance::calculateNortonEquivalentResistance(double time)
{//����֧·��ŵ�ٵ�Ч����
	nortonEquivalentResistance = 2*inductanceValue/deltaT;
}