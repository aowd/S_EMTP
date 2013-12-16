#ifndef PICTRLCOMP_H
#define PICTRLCOMP_H

#include "CtrlComponent.h"

class PICtrlComp : public CtrlComponent
{
public:
	PICtrlComp(int id,int inNode,int outNode,double pParam,double iParam);
	~PICtrlComp(){};

	virtual void initializeCtrlBranch();//初始化支路输入输出信号值
	virtual void calculateCtrlEquivalentParameter();//计算元件参数
	virtual void saveInNodeValue(double* ctrlNodeValue);//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual void calculateInitOutputValue(double time);
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);

	// PI控制器初始化
	virtual void initializePICtrl(double* PICtrlInitialValue, int& ptr);

	// PWM变流器平均模型相关
	virtual void storeInternalVariables(); // 存储内部变量，以便在迭代校正时恢复
	virtual void restoreInternalVariables(); // 迭代校正时恢复内部变量

	virtual void storeInternalVariables_pre(); // 预测前存储
	virtual void restoreInternalVariables_pre(); // 预测后恢复

public:
	int inNode;//输入节点编号
	int outNode;
	double inNodeValue;//输入信号值
	double outNodeValue;//输出信号值
	double kp;//比例系数
	double ti;//积分系数  x=kp(1+1/(ti*s))u
	double c,d,c1,d1;//模块参数
	double histCurrent;
	double histCurrent_1; // 平均模型迭代校正使用
	double histCurrent_2; // 平均模型预测使用
};

#endif