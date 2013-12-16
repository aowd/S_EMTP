#ifndef PULSEGENCOMP_H
#define PULSEGENCOMP_H

#include "CtrlComponent.h"

class PulseGenComp : public CtrlComponent
{
public:
	PulseGenComp(int id,int outNode,double freq,double maxValue,double minValue,double duty,double phase);
	~PulseGenComp(){};

	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual int checkCalCondition(int* nodeCalMark){return 1;};
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int outNode;
	double outNodeValue;
	double frequency;
	double maxValue;
	double minValue;
	double duty;
	double initPhase;
};

#endif