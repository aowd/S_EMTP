#ifndef CONTROLLEDVOLTAGESOURCE_H
#define CONTROLLEDVOLTAGESOURCE_H

#include "SimpleBranch.h"

class ControlledVoltageSource:public SimpleBranch {
public:
	ControlledVoltageSource(int id,int formNode,int toNode,double innerResistance,int CtrlSystemNode);
	~ControlledVoltageSource(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	virtual void setControlledVariable(double* ctrlNodeValue);//根据控制信号设定等值电压
	virtual void setControlledVariableForSwitch(double* ctrlNodeValue);//根据控制信号设定等值电压

private:
	double voltageValue;//电压值
	double innerResistance;//电源内阻，单位为欧
	int CtrlSystemNode;//对应的控制信号节点编号
};

#endif