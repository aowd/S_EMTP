#include "NodeVoltageMsrComp.h"

NodeVoltageMsrComp::NodeVoltageMsrComp(int id, int msrNode, int ctrlNode)
{
	type = 1;
	this->id = id;
	this->msrNode = msrNode;
	this ->ctrlNode = ctrlNode;
}

void NodeVoltageMsrComp::initializeBranch()
{
	measurands=0;
}

double NodeVoltageMsrComp::getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches)
{
	measurands = nodeVoltageVec(msrNode);//ȡ�ò����ڵ��ѹֵ
	return measurands;
}
