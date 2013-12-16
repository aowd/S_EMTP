#ifndef IMPULSEGENCOMP_H
#define IMPULSEGENCOMP_H

#include "CtrlComponent.h"

class ImpulseGenComp : public CtrlComponent
{
public:
	ImpulseGenComp(int id,int outNode,double freq,double firstTime,double specValue,double normValue);
	~ImpulseGenComp(){};

	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark){return 1;};
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int outNode;
	double outNodeValue;
	double frequency;
	double specValue;
	double normValue;
	double firstTime;
};

#endif