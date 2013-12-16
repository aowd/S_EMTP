#ifndef D2TTRANS_H
#define D2TTRANS_H

#include "CtrlComponent.h"

class D2TTrans : public CtrlComponent
{
public:
	D2TTrans(int id,int nodeD,int nodeQ,int nodeA,int nodeB,int nodeC);
	~D2TTrans(){};

	virtual void initializeCtrlBranch();//初始化支路输入输出信号值
	virtual void saveInNodeValue(double* ctrlNodeValue);//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[2];//输入节点编号：d,q
	int outNode[3];//输出节点编号： a,b,c
	double inNodeValue[2];//输入信号值
	double outNodeValue[3];//输出信号值
};

#endif