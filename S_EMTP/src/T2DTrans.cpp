#include "T2DTrans.h"

T2DTrans::T2DTrans(int id,int nodeA,int nodeB,int nodeC,int nodeD,int nodeQ)
{
	type =8;
	nPort=5;
	nInPort=3;
	nOutPort=2;
	this->id = id;
	inNode[0] = nodeA;
	inNode[1] = nodeB;
	inNode[2] = nodeC;
	outNode[0] = nodeD;
	outNode[1] = nodeQ;
}


void T2DTrans::initializeCtrlBranch()
{
	for (int i=0;i<3;i++)
	{
		inNodeValue[i] = 0;
	}
	for (int i=0;i<2;i++)
	{
		outNodeValue[i] = 0;
	}
}

void T2DTrans::saveInNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<3;i++)
	{
		inNodeValue[i]=ctrlNodeValue[inNode[i]-1];
	}
}

void T2DTrans::saveOutNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<2;i++)
	{
		ctrlNodeValue[outNode[i]-1]=outNodeValue[i];
	}
}

void T2DTrans::calculateOutputValue(double time)
{
	outNodeValue[0] = 2.0/3.0*(inNodeValue[0]-0.5*inNodeValue[1]-0.5*inNodeValue[2]);
	outNodeValue[1] = sqrt(1.0/3)*(inNodeValue[1]-inNodeValue[2]);
}

int T2DTrans::checkCalCondition(int* nodeCalMark)
{
	int calMark=1;
	for (int i=0;i<3;i++)
	{
		calMark = calMark & nodeCalMark[inNode[i]-1];
	}
	return calMark;
}

void T2DTrans::markOutputNode(int* nodeCalMark)
{
	for (int i=0;i<2;i++)
	{
		nodeCalMark[outNode[i]-1] = 1;
	}
}

void T2DTrans::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
