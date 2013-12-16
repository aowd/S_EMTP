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
		measurands = ((WoundInducGenerator*)tempCom)->wr;//取得绕线转子异步电机的转速
		return measurands;
	}
	else 	if (tempCom->type==10)
	{
		measurands = ((InducGenerator*)tempCom)->wr;//取得绕线转子异步电机的转速
		return measurands;
	}
	else
	{
		cout<<"对应支路不是异步电机元件"<<endl;
		return 0;
	}
}
