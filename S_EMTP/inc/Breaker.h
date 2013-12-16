#ifndef BREAKER_H
#define BREAKER_H

#include "SimpleBranch.h"

class Breaker:public SimpleBranch {
public:
	Breaker(){};//���캯��
	~Breaker(){};//��������

	virtual void calculateNortonEquivalentCurrent(double time){};//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time){};//����֧·��ŵ�ٵ�Ч����
	virtual bool checkSwitch(double time){return 0;};

	virtual bool getState(){return state;};//���ؿ���״̬
	virtual void switchIt();//�任����״̬
	virtual double getSwitchRatio(){return 0;};//���ؿ��ض�����Ĳ�ֵ��
	virtual void setControlledVariableForSwitch(TVectorD ctrlStateNodeValue, int* stateNode,int nStateNode,TVectorD ctrlInputNodeValue, int* inputNode,int nInputNode){};
protected:
	double onValue;//��ͨ����ֵ����λΪŷ
	double offValue;//�ضϵ���ֵ����λΪŷ
	int state;//����״̬
};

#endif