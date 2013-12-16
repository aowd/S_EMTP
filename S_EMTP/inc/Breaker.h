#ifndef BREAKER_H
#define BREAKER_H

#include "SimpleBranch.h"

class Breaker:public SimpleBranch {
public:
	Breaker(){};//构造函数
	~Breaker(){};//析构函数

	virtual void calculateNortonEquivalentCurrent(double time){};//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time){};//计算支路的诺顿等效电阻
	virtual bool checkSwitch(double time){return 0;};

	virtual bool getState(){return state;};//返回开关状态
	virtual void switchIt();//变换开关状态
	virtual double getSwitchRatio(){return 0;};//返回开关动作点的插值比
	virtual void setControlledVariableForSwitch(TVectorD ctrlStateNodeValue, int* stateNode,int nStateNode,TVectorD ctrlInputNodeValue, int* inputNode,int nInputNode){};
protected:
	double onValue;//导通电阻值，单位为欧
	double offValue;//关断电阻值，单位为欧
	int state;//开关状态
};

#endif