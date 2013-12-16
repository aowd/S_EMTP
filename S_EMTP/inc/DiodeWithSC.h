#ifndef DIODEWITHSC_H
#define DIODEWITHSC_H

#include "Breaker.h"


class DiodeWithSC:public Breaker{
public:
	DiodeWithSC(int id,int fromNode,int toNode,double onValue,double offValue,double Rd,double Cd,int state);
	~DiodeWithSC(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	virtual bool checkSwitch(double time);//检测开关动作

private:
	double forwardBreakoverVoltage;//正向击穿电压
	double reverseWithstandVoltage;//反向击穿电压
	double forwardVoltageDrop;//导通压降

	double Rd;//缓冲吸收电路中的电阻
	double Cd;//缓冲吸收电路中的电容
};

#endif