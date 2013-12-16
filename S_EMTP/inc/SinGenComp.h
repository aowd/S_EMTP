#ifndef SINGENCOMP_H
#define SINGENCOMP_H

#include "CtrlComponent.h"

class SinGenComp : public CtrlComponent
{
public:
	SinGenComp(int id,int outNode,double mag,double freq,double phase);
	~SinGenComp(){};

	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark){return 1;};
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int outNode;
	double outNodeValue;
	double frequency;
	double magnitude;
	double initPhase;
};

#endif