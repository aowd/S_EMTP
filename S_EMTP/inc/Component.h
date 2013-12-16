#ifndef COMPONENT_H
#define COMPONENT_H
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

class Component
{
public:
	Component(){isUserDef = 0;need_NEC=0;};
	virtual ~Component(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time){};//初始化支路电压电流
	virtual void setDetalT(double deltaT){this->deltaT = deltaT;};
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray){};//从节点电压数组中读支路两节点电压
	virtual void calculateBranchVoltage(){};//计算支路电压
	virtual void calculateBranchCurrent(){};//计算支路电流
	virtual void calculateNortonEquivalentCurrent(double time){};//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time){};//计算支路的诺顿等效电阻
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray){};//形成节点诺顿等效电流向量
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix){};//形成节点导纳阵
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter){};//保存支路电流
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter){};//保存支路电流
	virtual void saveMachineWr(double** machineWrMatrix, int& ptr, int counter){};//保存电机转速
	virtual double getNortonEquivalentResistance(){return 0;};
	virtual double getNortonEquivalentCurrent(){return 0;};
	virtual double getNortonEquivalentCurrent_1(){return 0;};
	virtual double getNortonEquivalentCurrent_2(){return 0;};

	virtual void interpolate(double ratio){};//给定比值，对支路的电压电流进行插值
	virtual void updateResult(int updateTypeNum){};//更新开关处理过程中用于存储结果的变量

	virtual bool checkSwitch(double time){return 0;};//检测开关是否需要动作
	virtual bool checkSwitch(int counter,TMatrixD &conductanceMatrix){return 0;};//for Inverter
	virtual bool getState(){return 0;};//返回开关状态
	virtual double getSwitchRatio(){return 0;};//返回开关动作点的插值比
	virtual void switchIt(){};//变换开关状态
	virtual void modifyConductanceMatrix(TMatrixD &conductanceMatrix){};//修正节点导纳阵
	virtual int getSwitchMode(){return 0;};

	virtual void setControlledVariable(double* ctrlNodeValue){};
	virtual void setControlledVariableForSwitch(double* ctrlNodeValue){};

	// PWM变流器专用，具体说明参见PWMConverter.h,xuyin,20121208
	virtual void initializeSwitchingInstants(double** PWMInitialMatrix, int nStep, int& ptr){}; // 初始化开关时刻
	virtual void ctrlNodeNumber_PWM(int* ctrlNodeNumber, int k){}; // 记录PWM变流对应的控制节点编号
	virtual void predictSwithcingInstants(){}; // 预测开关时刻
	virtual int predictSwithcingInstants(double** PWMPredictMatrix, int nStep, int& ptr){return 0;}; // 预测开关时刻
	virtual int correctSwithcingInstants(double** PWMCorrectMatrix, double tol, int nStep, int& ptr){return 0;}; // 校正开关时刻
	
	// 配合平均模型的迭代校正算法使用
	// 目前仅PWMConverter和SimpleBranch中添加了这两个函数，xuyin,20121208
	virtual void storeInternalVariables(){}; // 存储内部变量，以便在迭代校正时恢复
	virtual void restoreInternalVariables(){}; // 迭代校正时恢复内部变量

	// 电机初始化，xuyin, 20121226
	virtual void initializeGen(double** GenInitialMatrix, int& ptr){};
	virtual void getWindVelocityData(double** VwMatrix, int rows, int& ptr){};

	// xuyin, 20130224
	virtual void correctYne(double* tao_s){};

	//检测可变电阻状态变化
	virtual int timedResistanceStatusChange(double time){return 0;};

public:
	int id;
	string name;
	int type;
	int nPort;
	double deltaT;
	double branchCurrent;
	double branchCurrent_1;
	double branchCurrent_2;
	double branchVoltage;
	double branchVoltage_1;
	double branchVoltage_2;
	int *nodeNumber;

	//是否自定义元件
	bool isUserDef;

	//是否需要计算诺顿等值电流
	bool need_NEC;	

	// PWM变流器的开关动作时刻, xuyin, 20130224
	double tao_s[3];
};

#endif