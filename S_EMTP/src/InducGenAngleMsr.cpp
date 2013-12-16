#include "InducGenAngleMsr.h"

InducGenAngleMsr::InducGenAngleMsr(int id, int msrBranch, int ctrlNode)
{
	type = 5;
	this->id = id;
	this->msrBranch = msrBranch;
	this ->ctrlNode = ctrlNode;
}

void InducGenAngleMsr::initializeBranch()
{
	measurands=0;
}

double InducGenAngleMsr::getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches)
{
	Component* tempCom;
	tempCom = branches->at(msrBranch-1);
	if (tempCom->type==18)
	{
		measurands = ((WoundInducGenerator*)tempCom)->curAngle_s-((WoundInducGenerator*)tempCom)->curAngle_r;//ȡ������ת���첽�����ת��
		return measurands;
	}
	else 	if (tempCom->type==10)
	{
		measurands = ((InducGenerator*)tempCom)->curAngle_s - ((InducGenerator*)tempCom)->curAngle_r;//ȡ������ת���첽�����ת��
		return measurands;
	}
	else
	{
		cout<<"��Ӧ֧·�����첽���Ԫ��"<<endl;
		return 0;
	}
}
