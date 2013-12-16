#ifndef NODEVOLTAGEMSRCOMP_H
#define NODEVOLTAGEMSRCOMP_H

#include "MeasureComponent.h"

class NodeVoltageMsrComp : public MeasureComponent
{
public:
	NodeVoltageMsrComp(int id, int msrNode,int ctrlNode);
	~NodeVoltageMsrComp(){};

	virtual void initializeBranch();//��ʼ��֧·��ѹ����
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//��ò���ֵ

public:
	int msrNode;//�����ڵ���
};


#endif