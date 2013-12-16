#include "BranchVoltageMsrComp.h"

BranchVoltageMsrComp::BranchVoltageMsrComp(int id, int msrBranch, int ctrlNode)
{
	type = 2;
	this->id = id;
	this->msrBranch = msrBranch;
	this ->ctrlNode = ctrlNode;
}

void BranchVoltageMsrComp::initializeBranch()
{
	measurands=0;
}

double BranchVoltageMsrComp::getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches)
{
	Component* tempCom;
	tempCom = branches->at(msrBranch-1);
	measurands = tempCom->branchVoltage;//ȡ�ò���֧·��ѹֵ
	return measurands;
}
