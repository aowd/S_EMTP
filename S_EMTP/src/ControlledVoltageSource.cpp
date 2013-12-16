#include "ControlledVoltageSource.h"
#include <cmath>
#include <iostream>
using namespace std;

ControlledVoltageSource::ControlledVoltageSource(int id,int fromNode,int toNode,double innerResistance,int CtrlSystemNode)
{
	type=23;
	this->id=id;
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	this->innerResistance = innerResistance;
	this->CtrlSystemNode = CtrlSystemNode;
	nortonEquivalentResistance = innerResistance;
	need_NEC=1;
	isSeriesOrNot = isSeries();
}

void ControlledVoltageSource::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
	nortonEquivalentCurrent = voltageValue/innerResistance;
}

void ControlledVoltageSource::calculateNortonEquivalentResistance(double time)
{
	nortonEquivalentResistance = innerResistance;
}

void ControlledVoltageSource::setControlledVariable(double* ctrlNodeValue)
{
	voltageValue = ctrlNodeValue[CtrlSystemNode-1];
}

void ControlledVoltageSource::setControlledVariableForSwitch(double* ctrlNodeValue)
{
	voltageValue = ctrlNodeValue[CtrlSystemNode-1];
}