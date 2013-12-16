#ifndef INDUCGENANGLEMSR_H
#define INDUCGENANGLEMSR_H

#include "MeasureComponent.h"
#include "WoundInducGenerator.h"
#include "InducGenerator.h"

class InducGenAngleMsr : public MeasureComponent
{
public:
	InducGenAngleMsr(int id, int msrBranch,int ctrlNode);
	~InducGenAngleMsr(){};

	virtual void initializeBranch();//初始化测量输入信号
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//获取测量信号值

public:
	int msrBranch;//测量支路编号
};

#endif