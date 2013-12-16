#ifndef SAMPLER_H
#define SAMPLER_H

#include "CtrlComponent.h"

class Sampler : public CtrlComponent
{
public:
	Sampler(int id,int inNode,int pulseNode, int outNode);
	~Sampler(){};

	virtual void initializeCtrlBranch();//初始化支路输入输出信号值
	virtual void saveInNodeValue(double* ctrlNodeValue);//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual void calculateInitOutputValue(double time);
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);

	// PWM变流器平均模型相关
	virtual void storeInternalVariables(); // 存储内部变量，以便在迭代校正时恢复
	virtual void restoreInternalVariables(); // 迭代校正时恢复内部变量

	virtual void storeInternalVariables_pre(); // 预测前存储
	virtual void restoreInternalVariables_pre(); // 预测后恢复

public:
	int inNode;//输入节点编号
	int pulseNode;
	int outNode;
	double inNodeValue;//输入信号值
	double pulseValue;//采样控制信号值
	double outNodeValue;//输出信号值
	double histValue;
	double histValue_1; // 平均模型迭代校正使用
	double histValue_2; // 平均模型预测使用
};

#endif