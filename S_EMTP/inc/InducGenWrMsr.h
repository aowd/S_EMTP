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

	virtual void initializeBranch();//��ʼ�����������ź�
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//��ȡ�����ź�ֵ

public:
	int msrBranch;//����֧·���
};

#endif