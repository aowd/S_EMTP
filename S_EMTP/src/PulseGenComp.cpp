#include "PulseGenComp.h"

PulseGenComp::PulseGenComp(int id,int outNode,double freq,double maxValue,double minValue,double duty,double phase)
{
	type = 4;
	nPort=1;
	nInPort=0;
	nOutPort =1;
	this->id = id;
	this->outNode=outNode;
	this->frequency= freq;
	this->maxValue=maxValue;
	this->minValue=minValue;
	this->duty=duty;
	this->initPhase=phase;
}

void PulseGenComp::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1] = outNodeValue;
}

void PulseGenComp::calculateOutputValue(double time)
{
	double phase=time*frequency+initPhase/360;
	double rate = phase-(int)phase;
	if (rate<=duty)
	{
		outNodeValue = maxValue;
	} 
	else
	{
		outNodeValue =  minValue;
	}
}

void PulseGenComp::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void PulseGenComp::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
