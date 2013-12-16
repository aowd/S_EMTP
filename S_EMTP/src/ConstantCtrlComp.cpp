#include "ConstantCtrlComp.h"

ConstantCtrlComp::ConstantCtrlComp(int id,int outNode,double constant)
{
	type = 3;
	nPort=1;
	nInPort=0;
	nOutPort=1;
	this->id = id;
	this->outNode=outNode;
	this->outNodeValue=constant;
}


void ConstantCtrlComp::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1]=outNodeValue;
}

void ConstantCtrlComp::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void ConstantCtrlComp::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}