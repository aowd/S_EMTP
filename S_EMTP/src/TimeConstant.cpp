#include "TimeConstant.h"

TimeConstant::TimeConstant(int id,int outNode)
{
	type = 15;
	nPort=1;
	nInPort=0;
	nOutPort=1;
	this->id = id;
	this->outNode=outNode;
}

void TimeConstant::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1]=outNodeValue;
}

void TimeConstant::calculateOutputValue(double time)
{
	outNodeValue = time;
}

void TimeConstant::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void TimeConstant::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}