#ifndef INDUCTANCE_H
#define INDUCTANCE_H

#include "SimpleBranch.h"

class Inductance:public SimpleBranch {
public:
	Inductance(int id,int formNode,int toNode,double value);
	~Inductance(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����

private:
	double inductanceValue;//���ֵ����λΪ��
};

#endif