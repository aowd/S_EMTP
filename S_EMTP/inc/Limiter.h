#ifndef LIMITER_H
#define LIMITER_H

#include "CtrlComponent.h"

class Limiter : public CtrlComponent
{
public:
	Limiter(int id,int inNode,int outNode,double upLim,double downLim);
	~Limiter(){};

	virtual void initializeCtrlBranch();//��ʼ��֧·��������ź�ֵ
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode;//����ڵ���
	int outNode;
	double inNodeValue;//�����ź�ֵ
	double outNodeValue;//����ź�ֵ
	double upLim;//����ֵ
	double downLim;//����ֵ
};

#endif