#include "PCtrlComp.h"

PCtrlComp::PCtrlComp(int id,int inNode,int outNode,double pParam)
{
	type =7;
	nPort=2;
	nInPort=1;
	nOutPort=1;
	this->id = id;
	this->inNode = inNode;
	this->outNode = outNode;
	this->kp = pParam;
}


void PCtrlComp::initializeCtrlBranch()
{
	inNodeValue = 0;
	outNodeValue = 0;
}

void PCtrlComp::saveInNodeValue(double* ctrlNodeValue)
{
	inNodeValue=ctrlNodeValue[inNode-1];
}

void PCtrlComp::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1]=outNodeValue;
}

void PCtrlComp::calculateOutputValue(double time)
{
	outNodeValue = kp*inNodeValue;
}

int PCtrlComp::checkCalCondition(int* nodeCalMark)
{
	return nodeCalMark[inNode-1];
}

void PCtrlComp::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void PCtrlComp::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
