#ifndef TIMEDRESISTANCE_H
#define TIMEDRESISTANCE_H

#include "SimpleBranch.h"

class TimedResistance:public SimpleBranch {
public:
	TimedResistance(int id,int formNode,int toNode,double value, double Tstart, double Tend, double changeValue);
	~TimedResistance(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	virtual int timedResistanceStatusChange(double time);//判断电阻值是否发生改变

	virtual void storeInternalVariables(); // 存储内部变量，以便在迭代校正时恢复
	virtual void restoreInternalVariables(); // 迭代校正时恢复内部变量

private:
	double resistanceValue;//电阻值，单位为欧
	double nortonEquivalentResistance_bak;//保存诺顿等值电阻值，以便迭代时恢复
	double Tstart;//变化开始时间
	double Tend;//变化结束时间
	double changeValue;//变化电阻值
};

#endif