#ifndef S2RTRANS_H
#define S2RTRANS_H

#include "CtrlComponent.h"

class S2RTrans : public CtrlComponent
{
public:
	S2RTrans(int id,int nodeAlpha,int nodeBeta,int nodeRho,int nodeD,int nodeQ);
	~S2RTrans(){};

	virtual void initializeCtrlBranch();//初始化支路输入输出信号值
	virtual void saveInNodeValue(double* ctrlNodeValue);//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[3];//输入节点编号：alpha,beta,rho
	int outNode[2];//输出节点编号： d,q
	double inNodeValue[3];//输入信号值
	double outNodeValue[2];//输出信号值
};

#endif