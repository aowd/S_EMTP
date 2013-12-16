#ifndef CTRLCOMPONENT_H
#define CTRLCOMPONENT_H
#include <windows.h>
#ifdef WIN32
#include <w32pragma.h>
#endif
#include <Riostream.h>
#include <TMath.h>
#include <TMatrixDUtils.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompLU.h>

#include <string>

#define PI 3.141592653589793238462643383279

class CtrlComponent
{
public:
	CtrlComponent(){};
	virtual ~CtrlComponent(){};

	virtual void initializeCtrlBranch(){};//初始化控制支路输入输出信号
	virtual void calculateCtrlEquivalentParameter(){};//计算支路中各个参数
	virtual void saveInNodeValue(double* ctrlNodeValue){};//保存输入控制信号
	virtual void saveOutNodeValue(double* ctrlNodeValue){};//保存输出控制信号
	virtual void calculateOutputValue(double time){};//计算模块输出信号值
	virtual void calculateInitOutputValue(double time){};//计算模块输出信号值
	virtual int checkCalCondition(int* nodeCalMark){return 0;};
	virtual void markOutputNode(int* nodeCalMark){};
	virtual void setDeltaT(double deltaT){this->deltaT = deltaT;};

	// PWM变流器平均模型相关
	virtual void storeInternalVariables(){}; // 存储内部变量，以便在迭代校正时恢复
	virtual void restoreInternalVariables(){}; // 迭代校正时恢复内部变量

	virtual void storeInternalVariables_pre(){}; // 预测前存储
	virtual void restoreInternalVariables_pre(){}; // 预测后恢复

	// xuyin, 20121110, PI控制器初始化
	virtual void initializePICtrl(double* PICtrlInitialValue, int& ptr){};
	virtual void initializeDelay(double** DelayInitialValue, int& ptr, int nRows){};

public:
	int id;
	string name;
	int type;
	int nPort;//端口个数
	int nInPort;
	int nOutPort;
	double deltaT;
};

#endif