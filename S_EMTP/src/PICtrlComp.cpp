#include "PICtrlComp.h"

PICtrlComp::PICtrlComp(int id,int inNode,int outNode,double pParam, double iParam)
{
	type =2;
	nPort =2;
	nInPort = 1;
	nOutPort = 1;
	this->id = id;
	this->inNode = inNode;
	this->outNode = outNode;
	this->kp = pParam;
	this->ti = iParam;
}


void PICtrlComp::initializeCtrlBranch()
{
	inNodeValue = 0;
	outNodeValue = 0;
}

// PI��������ʼ��
void PICtrlComp::initializePICtrl(double* PICtrlInitialValue, int& ptr)
{
	inNodeValue = PICtrlInitialValue[ptr];
	outNodeValue = PICtrlInitialValue[ptr+1];
	ptr += 2;
}

void PICtrlComp::calculateCtrlEquivalentParameter()
{
	c = 2/deltaT*ti;
	c1=-2/deltaT*ti;
	/*d=1+2/deltaT*ti;
	d1=1-2/deltaT*ti;*/
	d=1+2/deltaT*ti*kp;
	d1=1-2/deltaT*ti*kp;
}


void PICtrlComp::saveInNodeValue(double* ctrlNodeValue)
{
	inNodeValue = ctrlNodeValue[inNode-1];
}

void PICtrlComp::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1] = outNodeValue;
}

void PICtrlComp::calculateOutputValue(double time)
{
	outNodeValue = (d*inNodeValue+histCurrent)/c;
	histCurrent = d1*inNodeValue - c1*outNodeValue;
}

void PICtrlComp::calculateInitOutputValue(double time)
{
	// outNodeValue = kp*inNodeValue;
	histCurrent = d1*inNodeValue - c1*outNodeValue;
}

int PICtrlComp::checkCalCondition(int* nodeCalMark)
{
	return nodeCalMark[inNode-1];
}

void PICtrlComp::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

// PWM������ƽ��ģ�����
// �洢�ڲ��������Ա��ڵ���У��ʱ�ָ���PI�������洢��ʷֵ����
void PICtrlComp::storeInternalVariables()
{
	histCurrent_1 = histCurrent;
}

// ����У��ʱ�ָ��ڲ�����
void PICtrlComp::restoreInternalVariables()
{
	histCurrent = histCurrent_1;
}

// Ԥ��ǰ�洢
void PICtrlComp::storeInternalVariables_pre()
{
	histCurrent_2 = histCurrent;
}

// Ԥ���ָ�
void PICtrlComp::restoreInternalVariables_pre()
{
	histCurrent = histCurrent_2;
} 