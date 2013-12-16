#ifndef S2RTRANS_H
#define S2RTRANS_H

#include "CtrlComponent.h"

class S2RTrans : public CtrlComponent
{
public:
	S2RTrans(int id,int nodeAlpha,int nodeBeta,int nodeRho,int nodeD,int nodeQ);
	~S2RTrans(){};

	virtual void initializeCtrlBranch();//��ʼ��֧·��������ź�ֵ
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[3];//����ڵ��ţ�alpha,beta,rho
	int outNode[2];//����ڵ��ţ� d,q
	double inNodeValue[3];//�����ź�ֵ
	double outNodeValue[2];//����ź�ֵ
};

#endif