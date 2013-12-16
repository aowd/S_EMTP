#ifndef PCTRLCOMP_H
#define PCTRLCOMP_H

#include "CtrlComponent.h"

class PCtrlComp : public CtrlComponent
{
public:
	PCtrlComp(int id,int inNode,int outNode,double pParam);
	~PCtrlComp(){};

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
	double kp;//����ϵ��
};

#endif