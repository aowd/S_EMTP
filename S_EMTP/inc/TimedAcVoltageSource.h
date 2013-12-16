#ifndef TIMEDACVOLTAGESOURCE_H
#define TIMEDACVOLTAGESOURCE_H

#include "SimpleBranch.h"

class TimedAcVoltageSource:public SimpleBranch {
public:
	TimedAcVoltageSource(int id,int fromNode,int toNode,double mag,double phase,double freq,double innerResistance,double Tstart, double Tend, double dropRatio);
	~TimedAcVoltageSource(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻

private:
	double magnitude;//电压幅值，单位为伏
	double initialAngle;//初相角，单位为度
	double frequency;//频率，单位为赫兹
	double innerResistance;//电源内阻，单位为欧
	double Tstart;//故障开始时间
	double Tend;//故障结束时间
	double dropRatio;//故障过程中电压与原幅值之比
};

#endif