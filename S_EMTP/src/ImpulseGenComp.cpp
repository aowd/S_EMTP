#include "ImpulseGenComp.h"

ImpulseGenComp::ImpulseGenComp(int id,int outNode,double freq,double firstTime, double specValue,double normValue)
{
	type = 16;
	nPort=1;
	nInPort=0;
	nOutPort=1;
	this->id = id;
	this->outNode=outNode;
	this->frequency= freq;
	this->specValue=specValue;
	this->normValue=normValue;
	this->firstTime=firstTime;
}

void ImpulseGenComp::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1]=outNodeValue;
}

void ImpulseGenComp::calculateOutputValue(double time)
{
	double phase = time*frequency;
	int phaseInt = (int) (phase+0.5);
	if (time<firstTime-deltaT/2)
	{
		outNodeValue = normValue;
	}
	else if (phase>phaseInt-frequency*deltaT/2 && phase<phaseInt+frequency*deltaT/2 )
	{
		outNodeValue = specValue;
	}
	else
	{
		outNodeValue = normValue;
	}
}

void ImpulseGenComp::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void ImpulseGenComp::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
