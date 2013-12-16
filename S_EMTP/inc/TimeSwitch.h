#ifndef TIMESWITCH_H
#define TIMESWITCH_H

#include "Breaker.h"

class TimeSwitch:public Breaker{
public:
	TimeSwitch(int id,int fromNode,int toNode,double onValue,double offValue,int state,int actionTimeCase);
	~TimeSwitch(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual bool checkSwitch(double time);//��⿪�ض���

private:
	TVectorD actionTimeArray;//���忪�ض���ʱ��
};

#endif