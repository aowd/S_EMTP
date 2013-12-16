#ifndef COMPARATOR_H
#define COMPARATOR_H

#include "CtrlComponent.h"

class Comparator : public CtrlComponent
{
public:
	Comparator(int id,int nodeA,int nodeB,int nodeC,int nodeD,double value1,double value2);
	~Comparator(){};

	virtual void initializeCtrlBranch();//初始化支路输入输出信号值
	virtual void saveInNodeValue(double* ctrlNodeValue);//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[2];//输入节点编号
	int outNode[2];//输出节点编号
	double inNodeValue[2];//输入信号值
	double outNodeValue[2];//输出信号值
	double value1;
	double value2;
};

#endif