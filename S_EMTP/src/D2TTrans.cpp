#include "D2TTrans.h"

D2TTrans::D2TTrans(int id,int nodeD,int nodeQ,int nodeA,int nodeB,int nodeC)
{
	type =9;
	nPort=5;
	nInPort=2;
	nOutPort=3;
	this->id = id;
	outNode[0] = nodeA;
	outNode[1] = nodeB;
	outNode[2] = nodeC;
	inNode[0] = nodeD;
	inNode[1] = nodeQ;
}


void D2TTrans::initializeCtrlBranch()
{
	for (int i=0;i<2;i++)
	{
		inNodeValue[i] = 0;
	}
	for (int i=0;i<3;i++)
	{
		outNodeValue[i] = 0;
	}
}

void D2TTrans::saveInNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<2;i++)
	{
		inNodeValue[i]=ctrlNodeValue[inNode[i]-1];
	}
}

void D2TTrans::saveOutNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<3;i++)
	{
		ctrlNodeValue[outNode[i]-1]=outNodeValue[i];
	}
}

void D2TTrans::calculateOutputValue(double time)
{
	outNodeValue[0] = inNodeValue[0];
	outNodeValue[1] = 0.5*( -inNodeValue[0] + sqrt(3.0)* inNodeValue[1]);
	outNodeValue[2] = 0.5*( -inNodeValue[0] - sqrt(3.0)* inNodeValue[1]);
}

int D2TTrans::checkCalCondition(int* nodeCalMark)
{
	int calMark=1;
	for (int i=0;i<2;i++)
	{
		calMark = calMark & nodeCalMark[inNode[i]-1];
	}
	return calMark;
}

void D2TTrans::markOutputNode(int* nodeCalMark)
{
	for (int i=0;i<3;i++)
	{
		nodeCalMark[outNode[i]-1] = 1;
	}
}

void D2TTrans::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
