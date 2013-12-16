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

// PWM变流器平均模型相关
// 存储内部变量，以便在迭代校正时恢复，PI控制器存储历史值即可
void Sampler::storeInternalVariables()
{
	histValue_1 = histValue;
}

// 迭代校正时恢复内部变量
void Sampler::restoreInternalVariables()
{
	histValue = histValue_1;
}

// 预测前存储
void Sampler::storeInternalVariables_pre()
{
	histValue_2 = histValue;
}

// 预测后恢复
void Sampler::restoreInternalVariables_pre()
{
	histValue = histValue_2;
} 