#ifndef TIMESWITCH_H
#define TIMESWITCH_H

#include "Breaker.h"

class TimeSwitch:public Breaker{
public:
	TimeSwitch(int id,int fromNode,int toNode,double onValue,double offValue,int state,int actionTimeCase);
	~TimeSwitch(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	virtual bool checkSwitch(double time);//检测开关动作

private:
	TVectorD actionTimeArray;//定义开关动作时刻
};

#endif