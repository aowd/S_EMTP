#include "R2STrans.h"

R2STrans::R2STrans(int id,int nodeD,int nodeQ,int nodeRho,int nodeAlpha,int nodeBeta)
{
	type =11;
	nPort=5;
	nInPort=3;
	nOutPort=2;
	this->id = id;
	inNode[0] = nodeD;
	inNode[1] = nodeQ;
	inNode[2] = nodeRho;
	outNode[0] = nodeAlpha;
	outNode[1] = nodeBeta;
}


void R2STrans::initializeCtrlBranch()
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

void R2STrans::saveInNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<3;i++)
	{
		inNodeValue[i]=ctrlNodeValue[inNode[i]-1];
	}
}

void R2STrans::saveOutNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<2;i++)
	{
		ctrlNodeValue[outNode[i]-1]=outNodeValue[i];
	}
}

void R2STrans::calculateOutputValue(double time)
{
	outNodeValue[0] = inNodeValue[0]*cos(inNodeValue[2])-inNodeValue[1]*sin(inNodeValue[2]);
	outNodeValue[1] = inNodeValue[0]*sin(inNodeValue[2])+inNodeValue[1]*cos(inNodeValue[2]);
}

int R2STrans::checkCalCondition(int* nodeCalMark)
{
	int calMark=1;
	for (int i=0;i<3;i++)
	{
		calMark = calMark & nodeCalMark[inNode[i]-1];
	}
	return calMark;
}

void R2STrans::markOutputNode(int* nodeCalMark)
{
	for (int i=0;i<2;i++)
	{
		nodeCalMark[outNode[i]-1] = 1;
	}
}

void R2STrans::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
