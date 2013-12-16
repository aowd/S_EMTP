#ifndef DELAY_H
#define DELAY_H

#include "CtrlComponent.h"

class Delay : public CtrlComponent
{
public:
	Delay(int id,int inNode,int outNode,double tDelay,int nSamples);
	~Delay(){};

	virtual void initializeCtrlBranch();//初始化支路输入输出信号值
	virtual void saveInNodeValue(double* ctrlNodeValue);//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue);//保存输出控制信号
	virtual void calculateOutputValue(double time);//计算模块输出信号值
	virtual void calculateInitOutputValue(double time);
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);

	//初始化函数
	virtual void initializeDelay(double** DelayInitialValue, int& ptr,int nRows);

	// PWM变流器平均模型相关
	virtual void storeInternalVariables(); // 存储内部变量，以便在迭代校正时恢复
	virtual void restoreInternalVariables(); // 迭代校正时恢复内部变量

	virtual void storeInternalVariables_pre(); // 预测前存储
	virtual void restoreInternalVariables_pre(); // 预测后恢复

public:
	int inNode;//输入节点编号
	int outNode;//输出节点编号
	double inNodeValue;//输入信号值
	double outNodeValue;//输出信号值
	double tDelay;//延时时间
	int nSamples;//最大采样值个数
	int positionPointer;//队列指针位置
	int positionPointer_1;
	int positionPointer_2;
	int nCount;//循环次数指针
	int nCount_1;
	int nCount_2;
	double* storeQueue;//保存数值队列
	double* storeQueue_1;
	double* storeQueue_2;
	double histValue;//保存历史输出值
	double histValue_1;
	double histValue_2;
};

#endif