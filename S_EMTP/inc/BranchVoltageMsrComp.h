#ifndef BRANCHVOLTAGEMSRCOMP_H
#define BRANCHVOLTAGEMSRCOMP_H

#include "MeasureComponent.h"

class BranchVoltageMsrComp : public MeasureComponent
{
public:
	BranchVoltageMsrComp(int id, int msrBranch, int ctrlNode);
	~BranchVoltageMsrComp(){};

	virtual void initializeBranch();//初始化支路电压电流
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//获得测量量

public:
	int msrBranch;//测量支路编号
};

#endif