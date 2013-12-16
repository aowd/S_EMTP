#ifndef TRIANGLEGENCOMP_H
#define TRIANGLEGENCOMP_H

#include "CtrlComponent.h"

class TriangleGenComp : public CtrlComponent
{
public:
	TriangleGenComp(int id,int outNode,double freq,double maxValue,double minValue,double duty,double phase);
	~TriangleGenComp(){};

	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark){return 1;};
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int outNode;
	double outNodeValue;
	double frequency;
	double maxValue;
	double minValue;
	double duty;
	double initPhase;
};

#endif