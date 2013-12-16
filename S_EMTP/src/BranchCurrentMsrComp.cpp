#include "BranchCurrentMsrComp.h"

BranchCurrentMsrComp::BranchCurrentMsrComp(int id, int msrBranch, int ctrlNode)
{
	type = 3;
	this->id = id;
	this->msrBranch = msrBranch;
	this ->ctrlNode = ctrlNode;
}

void BranchCurrentMsrComp::initializeBranch()
{
	measurands=0;
}

//double BranchCurrentMsrComp::getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches)
//{
//	Component* tempCom;
//	tempCom = branches->at(msrBranch-1);
//	measurands = tempCom->branchCurrent;//��ȡ����֧·����ֵ
//	return measurands;
//}

// xuyin, 20121222
double BranchCurrentMsrComp::getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches)
{
	Component* tempCom;
	tempCom = branches->at(msrBranch-1);
	measurands = tempCom->branchCurrent;//��ȡ����֧·����ֵ

	//measurands = branchCurrentArray[msrBranch-1]; //��ȡ����֧·����ֵ
	return measurands;
}
