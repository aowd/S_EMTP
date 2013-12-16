#ifndef T2DTRANS_H
#define T2DTRANS_H

#include "CtrlComponent.h"

class T2DTrans : public CtrlComponent
{
public:
	T2DTrans(int id,int nodeA,int nodeB,int nodeC,int nodeD,int nodeQ);
	~T2DTrans(){};

	virtual void initializeCtrlBranch();//初始化支路输入输出信号值
	virtual void saveInNodeValue(double* ctrlNodeValue);//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[3];//输入节点编号：a,b,c
	int outNode[2];//输出节点编号： d,q
	double inNodeValue[3];//输入信号值
	double outNodeValue[2];//输出信号值
};

#endif