#ifndef BATTERY_H
#define BATTERY_H

#include "SimpleBranch.h"

class Battery:public SimpleBranch {
public:
	Battery(int id,int formNode,int toNode,double initSOC);
	~Battery(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual void calculateParameters();//����ͨ�õ�Чģ���е����ϵ��
	virtual void calculateBranchCurrent();//����Ԫ������
	virtual void storeInternalVariables();
	virtual void restoreInternalVariables();

private:
	double innerResistance;//����ֵ
	double Qfull,Efull,Eexp,Qexp,Enom,Qnom,Inom;//��ȫ���㡢ָ����ĩ�˵㡢�����ĩ�˵���ز���
	double A,B,K,E0;//ͨ�õ�Чģ������ز���
	double initSOC;//��ʼSOCֵ
	double Qit,Qit_bak;//��������ֵ���䱸��
};

#endif