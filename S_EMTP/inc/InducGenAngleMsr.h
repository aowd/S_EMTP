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

	virtual void initializeBranch();//��ʼ�����������ź�
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//��ȡ�����ź�ֵ

public:
	int msrBranch;//����֧·���
};

#endif