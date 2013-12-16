#ifndef DIODE_H
#define DIODE_H

#include "Breaker.h"


class Diode:public Breaker{
public:
	Diode(int id,int fromNode,int toNode,double onValue,double offValue,int state);
	~Diode(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	virtual bool checkSwitch(double time);//检测开关动作
	virtual double getSwitchRatio();//返回开关动作点的插值比
	virtual void switchIt();//变换开关状态
	virtual void modifyConductanceMatrix(TMatrixD &conductanceMatrix);//修正节点导纳阵
	virtual int getSwitchMode(){return 0;};


private:
	double forwardBreakoverVoltage;//正向击穿电压
	double reverseWithstandVoltage;//反向击穿电压
	double forwardVoltageDrop;//导通压降
};

#endif