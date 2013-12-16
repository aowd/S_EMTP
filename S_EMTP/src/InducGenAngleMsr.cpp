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
		measurands = ((WoundInducGenerator*)tempCom)->curAngle_s-((WoundInducGenerator*)tempCom)->curAngle_r;//取得绕线转子异步电机的转速
		return measurands;
	}
	else 	if (tempCom->type==10)
	{
		measurands = ((InducGenerator*)tempCom)->curAngle_s - ((InducGenerator*)tempCom)->curAngle_r;//取得绕线转子异步电机的转速
		return measurands;
	}
	else
	{
		cout<<"对应支路不是异步电机元件"<<endl;
		return 0;
	}
}
