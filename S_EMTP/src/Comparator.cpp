#include "Comparator.h"

Comparator::Comparator(int id,int nodeA,int nodeB,int nodeC,int nodeD,double value1, double value2)
{
	type =13;
	nPort=4;
	nInPort=2;
	nOutPort=2;
	this->id = id;
	inNode[0] = nodeA;
	inNode[1] = nodeB;
	outNode[0] = nodeC;
	outNode[1] = nodeD;
	this->value1 = value1;
	this->value2 = value2;
}


void Comparator::initializeCtrlBranch()
{
	for (int i=0;i<2;i++)
	{
		inNodeValue[i] = 0;
	}
	for (int i=0;i<2;i++)
	{
		outNodeValue[i] = 0;
	}
}

void Comparator::saveInNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<2;i++)
	{
		inNodeValue[i]=ctrlNodeValue[inNode[i]-1];
	}
}

void Comparator::saveOutNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<2;i++)
	{
		ctrlNodeValue[outNode[i]-1]=outNodeValue[i];
	}
}

void Comparator::calculateOutputValue(double time)
{
	if (inNodeValue[0]>inNodeValue[1])
	{
		outNodeValue[0] = value1;
		outNodeValue[1] = value2;
	} 
	else
	{
		outNodeValue[0] = value2;
		outNodeValue[1] = value1;
	}
}

int Comparator::checkCalCondition(int* nodeCalMark)
{
	int calMark=1;
	for (int i=0;i<2;i++)
	{
		calMark = calMark & nodeCalMark[inNode[i]-1];
	}
	return calMark;
}

void Comparator::markOutputNode(int* nodeCalMark)
{
	for (int i=0;i<2;i++)
	{
		nodeCalMark[outNode[i]-1] = 1;
	}
}

void Comparator::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
