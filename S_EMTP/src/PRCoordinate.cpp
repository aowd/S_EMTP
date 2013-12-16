#include "PRCoordinate.h"

PRCoordinate::PRCoordinate(int id,int nodeA,int nodeB,int nodeC,int nodeD,int option)
{
	type =12;
	nPort=4;
	nInPort=2;
	nOutPort=2;
	this->id = id;
	inNode[0] = nodeA;
	inNode[1] = nodeB;
	outNode[0] = nodeC;
	outNode[1] = nodeD;
	this->option = option;
}


void PRCoordinate::initializeCtrlBranch()
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

void PRCoordinate::saveInNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<2;i++)
	{
		inNodeValue[i]=ctrlNodeValue[inNode[i]-1];
	}
}

void PRCoordinate::saveOutNodeValue(double* ctrlNodeValue)
{
	for (int i=0;i<2;i++)
	{
		ctrlNodeValue[outNode[i]-1]=outNodeValue[i];
	}
}

void PRCoordinate::calculateOutputValue(double time)
{
	if (option==1)//直角坐标->极坐标
	{
		outNodeValue[0] = sqrt(inNodeValue[0]*inNodeValue[0]+inNodeValue[1]*inNodeValue[1]);
		if (inNodeValue[0]==0)
		{
			if (inNodeValue[1]>0)
			{
				outNodeValue[1] = PI/2;
			} 
			else if (inNodeValue[1]==0)
			{
				outNodeValue[1] = 0;
			}
			else if (inNodeValue[1]<0)
			{
				outNodeValue[1] = -PI/2;
			}
		} 
		else
		{
			outNodeValue[1] = atan(inNodeValue[1]/inNodeValue[0]);
			if (inNodeValue[0]<0 && inNodeValue[1]<0)
			{
				outNodeValue[1] -=PI;
			} 
			else if (inNodeValue[0]<0)
			{
				outNodeValue[1] +=PI;
			}
		}
	} 
	else//极坐标->直角坐标
	{
		outNodeValue[0] = inNodeValue[0]*cos(inNodeValue[1]);
		outNodeValue[1] = inNodeValue[0]*sin(inNodeValue[1]);
	}
}

int PRCoordinate::checkCalCondition(int* nodeCalMark)
{
	int calMark=1;
	for (int i=0;i<2;i++)
	{
		calMark = calMark & nodeCalMark[inNode[i]-1];
	}
	return calMark;
}

void PRCoordinate::markOutputNode(int* nodeCalMark)
{
	for (int i=0;i<2;i++)
	{
		nodeCalMark[outNode[i]-1] = 1;
	}
}

void PRCoordinate::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
