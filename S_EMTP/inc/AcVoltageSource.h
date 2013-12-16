#ifndef ACVOLTAGESOURCE_H
#define ACVOLTAGESOURCE_H

#include "SimpleBranch.h"

class AcVoltageSource:public SimpleBranch {
public:
	AcVoltageSource(int id,int formNode,int toNode,double mag,double phase,double freq,double innerResistance);
	~AcVoltageSource(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����

private:
	double magnitude;//��ѹ��ֵ����λΪ��
	double initialAngle;//����ǣ���λΪ��
	double frequency;//Ƶ�ʣ���λΪ����
	double innerResistance;//��Դ���裬��λΪŷ
};

#endif