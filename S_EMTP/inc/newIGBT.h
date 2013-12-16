#ifndef NEWIGBT_H
#define NEWIGBT_H

#include "Breaker.h"


class newIGBT:public Breaker{
public:
	newIGBT(int id,int fromNode,int toNode,double onValue,double offValue,int CtrlSystemNode);
	~newIGBT(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	virtual bool checkSwitch(double time);//检测开关动作
	//virtual int setControlSignal(double time);
	//virtual void readControlSignal();
	virtual double getSwitchRatio();//返回开关动作点的插值比
	virtual void switchIt();//变换开关状态
	virtual void modifyConductanceMatrix(TMatrixD &conductanceMatrix);//修正节点导纳阵
	virtual int getSwitchMode();
	virtual void setControlledVariable(double* ctrlNodeValue);
	virtual void setControlledVariableForSwitch(double* ctrlNodeValue);//根据控制信号设定等值电压

private:
	double forwardBreakoverVoltage;//正向击穿电压
	double reverseWithstandVoltage;//反向击穿电压
	double forwardVoltageDrop;//导通压降
	int controlSignal;
	int controlSignal_1;
	int controlSignal_2;
	int controlChangeOrNot;
	int switchMode;
	//TMatrixD controlSignalMatrix;
	int CtrlSystemNode;
};

#endif