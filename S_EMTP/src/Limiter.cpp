#include "Limiter.h"

Limiter::Limiter(int id,int inNode,int outNode,double upLim, double downLim)
{
	type =14;
	nPort=2;
	nInPort=1;
	nOutPort=1;
	this->id = id;
	this->inNode = inNode;
	this->outNode = outNode;
	this->upLim = upLim;
	this->downLim = downLim;
}


void Limiter::initializeCtrlBranch()
{
	inNodeValue = 0;
	outNodeValue = 0;
}

void Limiter::saveInNodeValue(double* ctrlNodeValue)
{
	inNodeValue=ctrlNodeValue[inNode-1];
}

void Limiter::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1]=outNodeValue;
}

void Limiter::calculateOutputValue(double time)
{
	if (inNodeValue<downLim)
	{
		outNodeValue = downLim;
	} 
	else if (inNodeValue>upLim)
	{
		outNodeValue = upLim;
	}
	else
	{
		outNodeValue = inNodeValue;
	}
}

int Limiter::checkCalCondition(int* nodeCalMark)
{
	return nodeCalMark[inNode-1];
}

void Limiter::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

void Limiter::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}
