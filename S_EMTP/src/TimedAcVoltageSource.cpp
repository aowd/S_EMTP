#include "TimedAcVoltageSource.h"
#include <cmath>
#include <iostream>
using namespace std;

TimedAcVoltageSource::TimedAcVoltageSource(int id,int fromNode,int toNode,double mag,double phase,double freq,double innerResistance, double Tstart,double Tend,double dropRatio)
{
	type=20;
	this->id=id;
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	magnitude = mag;
	initialAngle = phase;
	frequency = freq;
	this->innerResistance = innerResistance;
	this->Tstart = Tstart;
	this->Tend = Tend;
	this-> dropRatio = dropRatio;
	need_NEC=1;
	isSeriesOrNot = isSeries();
}

void TimedAcVoltageSource::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
	if (time<Tstart-deltaT/2||time>Tend+deltaT/2)
	{
		nortonEquivalentCurrent = magnitude*sin(2*PI*frequency*time+initialAngle*PI/180)/innerResistance;
	}
	else
	{
		nortonEquivalentCurrent = dropRatio*magnitude*sin(2*PI*frequency*time+initialAngle*PI/180)/innerResistance;
	}
}

void TimedAcVoltageSource::calculateNortonEquivalentResistance(double time)
{
	nortonEquivalentResistance = innerResistance;
}
