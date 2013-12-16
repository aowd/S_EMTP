#ifndef CAPACITANCE_H
#define CAPACITANCE_H

#include "SimpleBranch.h"

class Capacitance:public SimpleBranch {
public:
	Capacitance(int id,int formNode,int toNode,double value);
	~Capacitance(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻

private:
	double capacitanceValue;//电容值，单位为法
};

#endif