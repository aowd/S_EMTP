#ifndef R2STRANS_H
#define R2STRANS_H

#include "CtrlComponent.h"

class R2STrans : public CtrlComponent
{
public:
	R2STrans(int id,int nodeD,int nodeQ,int nodeRho,int nodeAlpha,int nodeBeta);
	~R2STrans(){};

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