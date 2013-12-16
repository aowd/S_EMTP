#ifndef TIMECONSTANT_H
#define TIMECONSTANT_H

#include "CtrlComponent.h"

class TimeConstant : public CtrlComponent
{
public:
	TimeConstant(int id,int outNode);
	~TimeConstant(){};

	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual int checkCalCondition(int* nodeCalMark){return 1;};
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int outNode;
	double outNodeValue;
};

#endif