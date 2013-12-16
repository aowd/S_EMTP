#include "Resistance.h"

#include <iostream>
using namespace std;

Resistance::Resistance(int id,int fromNode,int toNode,double value)
{
	type=1;
	this->id=id;
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	resistanceValue =value;
	nortonEquivalentResistance = resistanceValue;

	isSeriesOrNot = isSeries();
}

void Resistance::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent_1 = 0;
	nortonEquivalentCurrent = 0;
}

void Resistance::calculateNortonEquivalentResistance(double time)
{
	nortonEquivalentResistance = resistanceValue;
}