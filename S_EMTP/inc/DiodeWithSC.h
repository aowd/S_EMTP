#ifndef DIODEWITHSC_H
#define DIODEWITHSC_H

#include "Breaker.h"


class DiodeWithSC:public Breaker{
public:
	DiodeWithSC(int id,int fromNode,int toNode,double onValue,double offValue,double Rd,double Cd,int state);
	~DiodeWithSC(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual bool checkSwitch(double time);//��⿪�ض���

private:
	double forwardBreakoverVoltage;//���������ѹ
	double reverseWithstandVoltage;//���������ѹ
	double forwardVoltageDrop;//��ͨѹ��

	double Rd;//�������յ�·�еĵ���
	double Cd;//�������յ�·�еĵ���
};

#endif