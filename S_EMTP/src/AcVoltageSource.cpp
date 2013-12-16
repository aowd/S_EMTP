#include "AcVoltageSource.h"
#include <cmath>
#include <iostream>
using namespace std;

AcVoltageSource::AcVoltageSource(int id,int fromNode,int toNode,double mag,double phase,double freq,double innerResistance)
{
	type=4;
	this->id=id;
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	magnitude = mag;
	initialAngle = phase;
	frequency = freq;
	this->innerResistance = innerResistance;
	nortonEquivalentResistance = innerResistance;
	need_NEC=1;
	isSeriesOrNot = isSeries();
}

void AcVoltageSource::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
	nortonEquivalentCurrent = magnitude*sin(2*PI*frequency*time+initialAngle*PI/180)/innerResistance;
}

void AcVoltageSource::calculateNortonEquivalentResistance(double time)
{
	nortonEquivalentResistance = innerResistance;
}
