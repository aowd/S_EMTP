#include "Sampler.h"

Sampler::Sampler(int id,int inNode,int pulseNode, int outNode)
{
	type =17;
	nPort =3;
	nInPort = 2;
	nOutPort = 1;
	this->id = id;
	this->inNode = inNode;
	this->outNode = outNode;
	this->pulseNode = pulseNode;
}


void Sampler::initializeCtrlBranch()
{
	inNodeValue = 0;
	pulseValue = 0;
	outNodeValue = 0;
}

void Sampler::saveInNodeValue(double* ctrlNodeValue)
{
	inNodeValue = ctrlNodeValue[inNode-1];
	pulseValue = ctrlNodeValue[pulseNode-1];
}

void Sampler::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1] = outNodeValue;
}

void Sampler::calculateOutputValue(double time)
{
	if (pulseValue==1)
	{
		histValue = inNodeValue;
	}
	outNodeValue = histValue;
}

void Sampler::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}

int Sampler::checkCalCondition(int* nodeCalMark)
{
	return nodeCalMark[inNode-1]&&nodeCalMark[pulseNode-1];
}

void Sampler::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

// PWM������ƽ��ģ�����
// �洢�ڲ��������Ա��ڵ���У��ʱ�ָ���PI�������洢��ʷֵ����
void Sampler::storeInternalVariables()
{
	histValue_1 = histValue;
}

// ����У��ʱ�ָ��ڲ�����
void Sampler::restoreInternalVariables()
{
	histValue = histValue_1;
}

// Ԥ��ǰ�洢
void Sampler::storeInternalVariables_pre()
{
	histValue_2 = histValue;
}

// Ԥ���ָ�
void Sampler::restoreInternalVariables_pre()
{
	histValue = histValue_2;
} 