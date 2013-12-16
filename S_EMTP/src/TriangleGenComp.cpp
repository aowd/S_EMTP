#include "TriangleGenComp.h"

TriangleGenComp::TriangleGenComp(int id,int outNode,double freq,double maxValue,double minValue,double duty,double phase)
{
	type = 5;
	nPort=1;
	nInPort=0;
	nOutPort=1;
	this->id = id;
	this->outNode=outNode;
	this->frequency= freq;
	this->maxValue=maxValue;
	this->minValue=minValue;
	this->duty=duty;
	this->initPhase=phase;
}

void TriangleGenComp::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1]=outNodeValue;
}

void TriangleGenComp::calculateOutputValue(double time)
{
	double phase=time*frequency+initPhase/360;
	double rate = phase-(int)phase;
	if (rate<=duty)
	{
		outNodeValue = minValue+(maxValue-minValue)/duty*rate;
	} 
	else
	{
		outNodeValue = maxValue+(minValue-maxValue)/(1-duty)*(rate-duty);
	}
}

void TriangleGenComp::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void TriangleGenComp::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
