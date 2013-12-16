#ifndef SIGMACTRLCOMPT_H
#define SIGMACTRLCOMPT_H

#include "CtrlComponent.h"

class SigmaCtrlComp : public CtrlComponent
{
public:
	SigmaCtrlComp(int id,int inNode1,int inNode2,int outNode,int option);
	~SigmaCtrlComp(){};

	virtual void initializeCtrlBranch();//初始化控制支路
	virtual void saveInNodeValue(double* ctrlNodeValue);//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[2];//输入节点编号集合
	int outNode;
	double inNodeValue[2];//输入信号值集合
	double outNodeValue;//输出信号值
	int option;//选择项，1为加和，2为做差、3为乘积、4为做商
};

#endif