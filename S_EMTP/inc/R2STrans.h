#ifndef R2STRANS_H
#define R2STRANS_H

#include "CtrlComponent.h"

class R2STrans : public CtrlComponent
{
public:
	R2STrans(int id,int nodeD,int nodeQ,int nodeRho,int nodeAlpha,int nodeBeta);
	~R2STrans(){};

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