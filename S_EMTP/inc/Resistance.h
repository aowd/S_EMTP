#ifndef RESISTANCE_H
#define RESISTANCE_H

#include "SimpleBranch.h"

class Resistance:public SimpleBranch {
public:
	Resistance(int id,int formNode,int toNode,double value);
	~Resistance(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����

private:
	double resistanceValue;//����ֵ����λΪŷ
};

#endif