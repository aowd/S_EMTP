#ifndef IGBT_H
#define IGBT_H

#include "Breaker.h"

class IGBT:public Breaker{
public:
	IGBT(int id,int fromNode,int toNode,double onValue,double offValue,int state,int actionTimeCase);
	~IGBT(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	virtual bool checkSwitch(double time);//检测开关动作

public:
	bool switchController(double time,int SwitchNnumber);
	bool camparePWM(double time,int SwitchNnumber);

public:
	double frequencyCW;
	double VmaxCW;
	double VoltageCW;
	double frequencyRef;
	double VmaxRef;
	double VoltageRef;
	int  IGBT_number;

};

#endif