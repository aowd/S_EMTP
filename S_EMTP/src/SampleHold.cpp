#include "SampleHold.h"

SampleHold::SampleHold(int id,int inNode, int selectNode, int outNode)
{
	type =20;
	nPort =3;
	nInPort = 2;
	nOutPort = 1;
	this->id = id;
	this->inNode = inNode;
	this->SelectNode = selectNode;
	this->outNode = outNode;
}

void SampleHold::initializeCtrlBranch()
{
	inNodeValue = 0;
	selectNodeValue =0;
	outNodeValue = 0;
	compStatus =0;
	compStatus_1 =0;
	compStatus_2 =0;
	histValue=0;
	histValue_1=0;
	histValue_2=0;
}

void SampleHold::saveInNodeValue(double* ctrlNodeValue)
{
	inNodeValue = ctrlNodeValue[inNode-1];
	selectNodeValue = ctrlNodeValue[SelectNode-1];
}

void SampleHold::saveOutNodeValue(double* ctrlNodeValue)
{
	ctrlNodeValue[outNode-1] = outNodeValue;
}

void SampleHold::calculateOutputValue(double time)
{
	if (selectNodeValue==0)
	{
		compStatus = 0;
		outNodeValue = inNodeValue;
	}
	else if (compStatus==0)
	{
		histValue = inNodeValue;
		compStatus =1;
		outNodeValue = histValue;
	}
	else
	{
		outNodeValue = histValue;
	}
}

void SampleHold::calculateInitOutputValue(double time)
{
	calculateOutputValue(time);
}

int SampleHold::checkCalCondition(int* nodeCalMark)
{
	return nodeCalMark[inNode-1]&&nodeCalMark[SelectNode-1];
}

void SampleHold::markOutputNode(int* nodeCalMark)
{
	nodeCalMark[outNode-1] = 1;
}

// PWM������ƽ��ģ�����
// �洢�ڲ��������Ա��ڵ���У��ʱ�ָ���PI�������洢��ʷֵ����
void SampleHold::storeInternalVariables()
{
	compStatus_1=compStatus;
	histValue_1= histValue;
}

// ����У��ʱ�ָ��ڲ�����
void SampleHold::restoreInternalVariables()
{
	compStatus = compStatus_1;
	histValue = histValue_1;
}

// Ԥ��ǰ�洢
void SampleHold::storeInternalVariables_pre()
{
	compStatus_2=compStatus;
	histValue_2 = histValue;
}

// Ԥ���ָ�
void SampleHold::restoreInternalVariables_pre()
{
	compStatus=compStatus_2;
	histValue = histValue_2;
} 