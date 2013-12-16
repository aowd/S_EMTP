#ifndef INDUCGENWRMSR_H
#define INDUCGENWRMSR_H

#include "MeasureComponent.h"
#include "WoundInducGenerator.h"
#include "InducGenerator.h"

class InducGenWrMsr : public MeasureComponent
{
public:
	InducGenWrMsr(int id, int msrBranch,int ctrlNode);
	~InducGenWrMsr(){};

	virtual void initializeBranch();//初始化测量输入信号
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//获取测量信号值

public:
	int msrBranch;//测量支路编号
};

#endif