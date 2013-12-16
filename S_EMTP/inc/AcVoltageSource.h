#ifndef ACVOLTAGESOURCE_H
#define ACVOLTAGESOURCE_H

#include "SimpleBranch.h"

class AcVoltageSource:public SimpleBranch {
public:
	AcVoltageSource(int id,int formNode,int toNode,double mag,double phase,double freq,double innerResistance);
	~AcVoltageSource(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻

private:
	double magnitude;//电压幅值，单位为伏
	double initialAngle;//初相角，单位为度
	double frequency;//频率，单位为赫兹
	double innerResistance;//电源内阻，单位为欧
};

#endif