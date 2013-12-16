#ifndef CAPACITANCE_H
#define CAPACITANCE_H

#include "SimpleBranch.h"

class Capacitance:public SimpleBranch {
public:
	Capacitance(int id,int formNode,int toNode,double value);
	~Capacitance(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����

private:
	double capacitanceValue;//����ֵ����λΪ��
};

#endif