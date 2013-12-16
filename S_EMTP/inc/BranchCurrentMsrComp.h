#ifndef BRANCHCURRENTMSRCOMP_H
#define BRANCHCURRENTMSRCOMP_H

#include "MeasureComponent.h"

class BranchCurrentMsrComp : public MeasureComponent
{
public:
	BranchCurrentMsrComp(int id, int msrBranch,int ctrlNode);
	~BranchCurrentMsrComp(){};

	virtual void initializeBranch();//初始化测量输入信号
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//获取测量信号值

public:
	int msrBranch;//测量支路编号
};

#endif