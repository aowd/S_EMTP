#ifndef RESISTANCE_H
#define RESISTANCE_H

#include "SimpleBranch.h"

class Resistance:public SimpleBranch {
public:
	Resistance(int id,int formNode,int toNode,double value);
	~Resistance(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻

private:
	double resistanceValue;//电阻值，单位为欧
};

#endif