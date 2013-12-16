#ifndef BRANCHCURRENTMSRCOMP_H
#define BRANCHCURRENTMSRCOMP_H

#include "MeasureComponent.h"

class BranchCurrentMsrComp : public MeasureComponent
{
public:
	BranchCurrentMsrComp(int id, int msrBranch,int ctrlNode);
	~BranchCurrentMsrComp(){};

	virtual void initializeBranch();//��ʼ�����������ź�
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//��ȡ�����ź�ֵ

public:
	int msrBranch;//����֧·���
};

#endif