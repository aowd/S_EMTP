#ifndef CONSTANTCTRLCOMP_H
#define CONSTANTCTRLCOMP_H

#include "CtrlComponent.h"

class ConstantCtrlComp : public CtrlComponent
{
public:
	ConstantCtrlComp(int id,int outNode,double constant);
	~ConstantCtrlComp(){};

	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual int checkCalCondition(int* nodeCalMark){return 1;};
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int outNode;
	double outNodeValue;
};

#endif