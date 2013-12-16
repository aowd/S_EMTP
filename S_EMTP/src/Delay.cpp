#include "Delay.h"

Delay::Delay(int id,int inNode,int outNode,double tDelay, int nSamples)
{
	type =18;
	nPort =2;
	nInPort = 1;
	nOutPort = 1;
	this->id = id;
	this->inNode = inNode;
	this->outNode = outNode;
	this->tDelay = tDelay;
	this->nSamples = nSamples;
}

void Delay::initializeCtrlBranch()
{
	inNodeValue = 0;
	outNodeValue = 0;
	positionPointer=0;
	positionPointer_1=0;
	positionPointer_2=0;
	int tempInt = (int)(tDelay/deltaT+0.5);
	nSamples = (tempInt>nSamples)?nSamples:tempInt;
	nCount=0;
	nCount_1=(int)(tDelay/deltaT+0.5)/nSamples-1;
	nCount_2=(int)(tDelay/deltaT+0.5)/nSamples-1;
	storeQueue = new double[nSamples];
	storeQueue_1 = new double[nSamples];
	storeQueue_2 = new double[nSamples];
	for (int i=0;i<nSamples;i++)
	{
		storeQueue[i]=0;
		storeQueue_1[i]=0;
		storeQueue_2[i]=0;
	}
	histValue=0;
	histValue_1=0;
	histValue_2=0;
}

void Delay::saveInNodeValue(double* ctrlNodeValue)
{
	inNodeValue = ctrlNodeValue[inNode-1];
}

void Delay::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1] = outNodeValue;
}

void Delay::calculateOutputValue(double time)
{
	nCount++;
	if (nCount==(int)(tDelay/deltaT+0.5)/nSamples)
	{
		nCount=0;
		histValue = storeQueue[positionPointer];
		storeQueue[positionPointer] = inNodeValue;	
		positionPointer++;
		positionPointer = (positionPointer==nSamples)?0:positionPointer;
	}
	outNodeValue = histValue;
}

void Delay::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}

void Delay::initializeDelay(double** DelayInitialValue, int& ptr, int nRows)
{
	nCount = 0;
	positionPointer= 0;
	for (int i=0;i<(int)(tDelay/deltaT+0.5);i++)
	{
		nCount++;
		if (nCount==(int)(tDelay/deltaT+0.5)/nSamples)
		{
			nCount=0;
			storeQueue[positionPointer] = DelayInitialValue[nRows-(int)(tDelay/deltaT+0.5)+i][ptr];	
			positionPointer++;
		}
	}
	ptr ++;
	nCount = 0;
	positionPointer = 0;
}

int Delay::checkCalCondition(int* nodeCalMark)
{
	return nodeCalMark[inNode-1];
}

void Delay::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

// PWM������ƽ��ģ�����
// �洢�ڲ��������Ա��ڵ���У��ʱ�ָ���PI�������洢��ʷֵ����
void Delay::storeInternalVariables()
{
	positionPointer_1 = positionPointer;
	nCount_1 = nCount;
	for (int i=0;i<nSamples;i++)
	{
		storeQueue_1[i] = storeQueue[i];
	}
	histValue_1= histValue;
}

// ����У��ʱ�ָ��ڲ�����
void Delay::restoreInternalVariables()
{
	positionPointer = positionPointer_1;
	nCount = nCount_1;
	for (int i=0;i<nSamples;i++)
	{
		storeQueue[i] = storeQueue_1[i];
	}
	histValue = histValue_1;
}

// Ԥ��ǰ�洢
void Delay::storeInternalVariables_pre()
{
	positionPointer_2 = positionPointer;
	nCount_2 = nCount;
	for (int i=0;i<nSamples;i++)
	{
		storeQueue_2[i] = storeQueue[i];
	}
	histValue_2 = histValue;
}

// Ԥ���ָ�
void Delay::restoreInternalVariables_pre()
{
	positionPointer = positionPointer_2;
	nCount = nCount_2;
	for (int i=0;i<nSamples;i++)
	{
		storeQueue[i] = storeQueue_2[i];
	}
	histValue = histValue_2;
} 