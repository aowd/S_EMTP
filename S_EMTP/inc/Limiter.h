#ifndef LIMITER_H
#define LIMITER_H

#include "CtrlComponent.h"

class Limiter : public CtrlComponent
{
public:
	Limiter(int id,int inNode,int outNode,double upLim,double downLim);
	~Limiter(){};

	virtual void initializeCtrlBranch();//初始化支路输入输出信号值
	virtual void saveInNodeValue(double* ctrlNodeValue);//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode;//输入节点编号
	int outNode;
	double inNodeValue;//输入信号值
	double outNodeValue;//输出信号值
	double upLim;//上限值
	double downLim;//下限值
};

#endif