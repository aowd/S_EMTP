#ifndef INDUCTANCE_H
#define INDUCTANCE_H

#include "SimpleBranch.h"

class Inductance:public SimpleBranch {
public:
	Inductance(int id,int formNode,int toNode,double value);
	~Inductance(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻

private:
	double inductanceValue;//电感值，单位为亨
};

#endif