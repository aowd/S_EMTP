#ifndef CONTROLLEDVOLTAGESOURCE_H
#define CONTROLLEDVOLTAGESOURCE_H

#include "SimpleBranch.h"

class ControlledVoltageSource:public SimpleBranch {
public:
	ControlledVoltageSource(int id,int formNode,int toNode,double innerResistance,int CtrlSystemNode);
	~ControlledVoltageSource(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual void setControlledVariable(double* ctrlNodeValue);//���ݿ����ź��趨��ֵ��ѹ
	virtual void setControlledVariableForSwitch(double* ctrlNodeValue);//���ݿ����ź��趨��ֵ��ѹ

private:
	double voltageValue;//��ѹֵ
	double innerResistance;//��Դ���裬��λΪŷ
	int CtrlSystemNode;//��Ӧ�Ŀ����źŽڵ���
};

#endif