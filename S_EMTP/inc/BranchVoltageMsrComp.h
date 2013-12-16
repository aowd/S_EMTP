#ifndef BRANCHVOLTAGEMSRCOMP_H
#define BRANCHVOLTAGEMSRCOMP_H

#include "MeasureComponent.h"

class BranchVoltageMsrComp : public MeasureComponent
{
public:
	BranchVoltageMsrComp(int id, int msrBranch, int ctrlNode);
	~BranchVoltageMsrComp(){};

	virtual void initializeBranch();//��ʼ��֧·��ѹ����
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//��ò�����

public:
	int msrBranch;//����֧·���
};

#endif