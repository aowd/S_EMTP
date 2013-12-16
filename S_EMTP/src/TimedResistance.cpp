#include "TimedResistance.h"

#include <iostream>
using namespace std;

TimedResistance::TimedResistance(int id,int fromNode,int toNode,double value, double Tstart, double Tend, double changeValue)
{
	type=27;
	this->id=id;
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	resistanceValue =value;
	this-> Tstart = Tstart;
	this-> Tend = Tend;
	this->changeValue = changeValue;
	need_NEC = 1;
	isSeriesOrNot = isSeries();
}

void TimedResistance::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent_1 = 0;
	nortonEquivalentCurrent = 0;
}

void TimedResistance::calculateNortonEquivalentResistance(double time)
{
	if (time<Tstart-deltaT/2||time>Tend-deltaT/2)
	{
		nortonEquivalentResistance = resistanceValue;
	}
	else
	{
		nortonEquivalentResistance = changeValue;
	}
}

int TimedResistance::timedResistanceStatusChange(double time)
{
	if ((time>=Tstart-deltaT/2&&time<Tstart+deltaT/2)||(time>=Tend-deltaT/2&&time<Tend+deltaT/2))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void TimedResistance::storeInternalVariables()
{
	nortonEquivalentResistance_bak=nortonEquivalentResistance;
}
void TimedResistance::restoreInternalVariables()
{
	nortonEquivalentResistance = nortonEquivalentResistance_bak;
}
