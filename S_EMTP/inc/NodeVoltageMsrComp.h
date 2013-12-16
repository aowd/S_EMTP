#ifndef NODEVOLTAGEMSRCOMP_H
#define NODEVOLTAGEMSRCOMP_H

#include "MeasureComponent.h"

class NodeVoltageMsrComp : public MeasureComponent
{
public:
	NodeVoltageMsrComp(int id, int msrNode,int ctrlNode);
	~NodeVoltageMsrComp(){};

	virtual void initializeBranch();//初始化支路电压电流
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches);//获得测量值

public:
	int msrNode;//测量节点编号
};


#endif