#include "Selector.h"

Selector::Selector(int id,int inNodeA,int inNodeB, int selectNode, int outNode,double selectAValue)
{
	type =19;
	nPort=4;
	nInPort=3;
	nOutPort=1;
	this->id = id;
	inNode = new int[nInPort];
	inNodeValue = new double [nInPort];
	inNode[0] = inNodeA;
	inNode[1] = inNodeB;
	inNode[2] = selectNode;
	this->outNode = outNode;
	this->selectAValue = selectAValue;
}


void Selector::initializeCtrlBranch()
{
	for (int i=0;i<nInPort;i++)
	{
		inNodeValue[i] = 0;
	}
	outNodeValue = 0;
}

void Selector::saveInNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<nInPort;i++)
	{
		inNodeValue[i] = ctrlNodeValue[inNode[i]-1];
	}
}

void Selector::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1]=outNodeValue;
}

void Selector::calculateOutputValue(double time)
{
	outNodeValue = (inNodeValue[2]==selectAValue)?inNodeValue[0]:inNodeValue[1];
}

int Selector::checkCalCondition(int* nodeCalMark)
{
	return nodeCalMark[inNode[0]-1]&&nodeCalMark[inNode[1]-1]&&nodeCalMark[inNode[2]-1];
}

void Selector::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void Selector::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
