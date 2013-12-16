#ifndef SELECTOR_H
#define SELECTOR_H

#include "CtrlComponent.h"

class Selector : public CtrlComponent
{
public:
	Selector(int id,int inNodeA, int inNodeB, int selectNode, int outNode,double selectAValue);
	~Selector(){};

	virtual void initializeCtrlBranch();//��ʼ��֧·��������ź�ֵ
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int* inNode;//����ڵ���
	int outNode;
	double* inNodeValue;//�����ź�ֵ
	double outNodeValue;//����ź�ֵ
	double selectAValue;//����ϵ��
};

#endif