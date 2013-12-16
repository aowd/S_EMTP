#include "SigmaCtrlComp.h"

SigmaCtrlComp::SigmaCtrlComp(int id,int inNode1,int inNode2,int outNode,int option)
{
	type = 1;
	nPort=3;
	nInPort=2;
	nOutPort=1;
	this->id = id;
	inNode[0] = inNode1;
	inNode[1] = inNode2;
	this->outNode = outNode;
	this->option = option;
}


void SigmaCtrlComp::initializeCtrlBranch()
{
	inNodeValue[0] = 0;
	inNodeValue[1] = 0;
	outNodeValue = 0;
}

void SigmaCtrlComp::saveInNodeValue(double* ctrlNodeValue)
{
	inNodeValue[0] = ctrlNodeValue[inNode[0]-1];
	inNodeValue[1] = ctrlNodeValue[inNode[1]-1];
}

void SigmaCtrlComp::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1] = outNodeValue;
}

void SigmaCtrlComp::calculateOutputValue(double time)
{
	switch(option)
	{
		case 1:
			outNodeValue = inNodeValue[0] + inNodeValue[1];
			break;
		case 2:
			outNodeValue = inNodeValue[0] - inNodeValue[1];
			break;
		case 3:
			outNodeValue = inNodeValue[0]*inNodeValue[1];
			break;
		case 4:
			outNodeValue = inNodeValue[0]/inNodeValue[1];
			break;
		default:
			cout<<"error in SigmaCtrlComp option";
	}
}

int SigmaCtrlComp::checkCalCondition(int* nodeCalMark)
{
	return nodeCalMark[inNode[0]-1] & nodeCalMark[inNode[1]-1];
}

void SigmaCtrlComp::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void SigmaCtrlComp::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}