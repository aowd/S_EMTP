#ifndef THREEPHASERECTIFIERBRIDGE_H
#define THREEPHASERECTIFIERBRIDGE_H

#include "Component.h"

class ThreePhaseRectifierBridge : public Component
{
public:

	ThreePhaseRectifierBridge(int id,int *ACNode,int *DCNode,int *state,double onValue,double offValue);
	~ThreePhaseRectifierBridge(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);
	virtual void calculateBranchVoltage();
	virtual void calculateBranchCurrent();
	virtual void calculateNortonEquivalentCurrent(double time);
	virtual void calculateNortonEquivalentResistance(double time);
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);
	virtual bool checkSwitch(double time);
	virtual void interpolate(double ratio);

private:
	void numberdiode();//给二极管编号
	void calculateDiodeVoltage();//计算各个二极管两端的电压
	void calculateDiodeCurrent();//计算各个二极管上通过的电流
	void calculatebranchCurrentArray();//计算四个支路上的电流
	void checkbridgemode();//判断整流桥状态：6+1（断续）
	void modifybridgemode();//修正整流桥状态
	void modifydiodemode();//修正二极管状态
	//int getState();
	//double getratio();

private:	
	double forwardBreakoverVoltage;//正向击穿电压
	double reverseWithstandVoltage;//反向击穿电压
	double forwardVoltageDrop;//导通压降
	double nortonEquivalentCurrent;
	double nodeVoltage[5];
	double diodeVoltage[6];
	//double diodeVoltage_1[6];
	double diodeCurrent[6];
	//double diodeCurrent_1[6];
	double branchCurrentArray[6];
	
	double nortonEquivalentResistance[6];
	double onValue;
	double offValue;
	double ratio_1;
	double switchtime;

	int diodenodeNo[6][2];
	int diodeState[6];
//	int nodeNumber[5];//这个是元件有的，尽量不要在头文件中定义变量大小
	int bridgeMode;
	int bridgeMode_tmp;
	int bridgeMode_1;
	int statechange[2];
	//int on_off;
};
#endif
