#ifndef SAMPLEHOLD_H
#define SAMPLEHOLD_H

#include "CtrlComponent.h"

class SampleHold : public CtrlComponent
{
public:
	SampleHold(int id,int inNode, int selectNode, int outNode);
	~SampleHold(){};

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
	int SelectNode;//选择节点编号
	int outNode;//输出节点编号
	double inNodeValue;//输入信号值
	double selectNodeValue;//采样信号值
	double outNodeValue;//输出信号值
	int compStatus;//0代表采样，1代表保持
	int compStatus_1;
	int compStatus_2;
	double histValue;//保存历史输出值
	double histValue_1;
	double histValue_2;
};

#endif