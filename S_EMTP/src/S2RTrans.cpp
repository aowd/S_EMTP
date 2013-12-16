#include "S2RTrans.h"

S2RTrans::S2RTrans(int id,int nodeAlpha,int nodeBeta,int nodeRho,int nodeD,int nodeQ)
{
	type =10;
	nPort=5;
	nInPort=3;
	nOutPort=2;
	this->id = id;
	inNode[0] = nodeAlpha;
	inNode[1] = nodeBeta;
	inNode[2] = nodeRho;
	outNode[0] = nodeD;
	outNode[1] = nodeQ;
}


void S2RTrans::initializeCtrlBranch()
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

void S2RTrans::saveInNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<3;i++)
	{
		inNodeValue[i]=ctrlNodeValue[inNode[i]-1];
	}
}

void S2RTrans::saveOutNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<2;i++)
	{
		ctrlNodeValue[outNode[i]-1]=outNodeValue[i];
	}
}

void S2RTrans::calculateOutputValue(double time)
{
	outNodeValue[0] = inNodeValue[0]*cos(inNodeValue[2])+inNodeValue[1]*sin(inNodeValue[2]);
	outNodeValue[1] = -inNodeValue[0]*sin(inNodeValue[2])+inNodeValue[1]*cos(inNodeValue[2]);
}

int S2RTrans::checkCalCondition(int* nodeCalMark)
{
	int calMark=1;
	for (int i=0;i<3;i++)
	{
		calMark = calMark & nodeCalMark[inNode[i]-1];
	}
	return calMark;
}

void S2RTrans::markOutputNode(int* nodeCalMark)
{
	for (int i=0;i<2;i++)
	{
		nodeCalMark[outNode[i]-1] = 1;
	}
}

void S2RTrans::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
