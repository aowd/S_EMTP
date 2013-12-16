#include "SinGenComp.h"

SinGenComp::SinGenComp(int id,int outNode,double mag,double freq,double phase)
{
	type = 6;
	nPort=1;
	nInPort=0;
	nOutPort=1;
	this->id = id;
	this->outNode=outNode;
	this->frequency= freq;
	this->magnitude=mag;
	this->initPhase=phase;
}

void SinGenComp::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1]=outNodeValue;
}

void SinGenComp::calculateOutputValue(double time)
{
	outNodeValue =magnitude*sin(2*PI*frequency*time+initPhase*PI/180);
}

void SinGenComp::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void SinGenComp::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}