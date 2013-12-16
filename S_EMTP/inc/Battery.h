#ifndef BATTERY_H
#define BATTERY_H

#include "SimpleBranch.h"

class Battery:public SimpleBranch {
public:
	Battery(int id,int formNode,int toNode,double initSOC);
	~Battery(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	virtual void calculateParameters();//计算通用等效模型中的相关系数
	virtual void calculateBranchCurrent();//计算元件电流
	virtual void storeInternalVariables();
	virtual void restoreInternalVariables();

private:
	double innerResistance;//内阻值
	double Qfull,Efull,Eexp,Qexp,Enom,Qnom,Inom;//完全充电点、指数区末端点、标称区末端点相关参数
	double A,B,K,E0;//通用等效模型中相关参数
	double initSOC;//初始SOC值
	double Qit,Qit_bak;//已用容量值及其备份
};

#endif