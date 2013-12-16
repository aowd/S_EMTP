#ifndef COMPARATOR_H
#define COMPARATOR_H

#include "CtrlComponent.h"

class Comparator : public CtrlComponent
{
public:
	Comparator(int id,int nodeA,int nodeB,int nodeC,int nodeD,double value1,double value2);
	~Comparator(){};

	virtual void initializeCtrlBranch();//��ʼ��֧·��������ź�ֵ
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[2];//����ڵ���
	int outNode[2];//����ڵ���
	double inNodeValue[2];//�����ź�ֵ
	double outNodeValue[2];//����ź�ֵ
	double value1;
	double value2;
};

#endif