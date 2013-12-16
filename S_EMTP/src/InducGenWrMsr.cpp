#include "InducGenWrMsr.h"

InducGenWrMsr::InducGenWrMsr(int id, int msrBranch, int ctrlNode)
{
	type = 4;
	this->id = id;
	this->msrBranch = msrBranch;
	this ->ctrlNode = ctrlNode;
}

void InducGenWrMsr::initializeBranch()
{
	measurands=0;
}

double InducGenWrMsr::getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches)
{
	Component* tempCom;
	tempCom = branches->at(msrBranch-1);
	if (tempCom->type==18)
	{
		measurands = ((WoundInducGenerator*)tempCom)->wr;//ȡ������ת���첽�����ת��
		return measurands;
	}
	else 	if (tempCom->type==10)
	{
		measurands = ((InducGenerator*)tempCom)->wr;//ȡ������ת���첽�����ת��
		return measurands;
	}
	else
	{
		cout<<"��Ӧ֧·�����첽���Ԫ��"<<endl;
		return 0;
	}
}
