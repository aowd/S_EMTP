#ifndef TIMEDRESISTANCE_H
#define TIMEDRESISTANCE_H

#include "SimpleBranch.h"

class TimedResistance:public SimpleBranch {
public:
	TimedResistance(int id,int formNode,int toNode,double value, double Tstart, double Tend, double changeValue);
	~TimedResistance(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual int timedResistanceStatusChange(double time);//�жϵ���ֵ�Ƿ����ı�

	virtual void storeInternalVariables(); // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(); // ����У��ʱ�ָ��ڲ�����

private:
	double resistanceValue;//����ֵ����λΪŷ
	double nortonEquivalentResistance_bak;//����ŵ�ٵ�ֵ����ֵ���Ա����ʱ�ָ�
	double Tstart;//�仯��ʼʱ��
	double Tend;//�仯����ʱ��
	double changeValue;//�仯����ֵ
};

#endif