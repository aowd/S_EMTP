#ifndef T2DTRANS_H
#define T2DTRANS_H

#include "CtrlComponent.h"

class T2DTrans : public CtrlComponent
{
public:
	T2DTrans(int id,int nodeA,int nodeB,int nodeC,int nodeD,int nodeQ);
	~T2DTrans(){};

	virtual void initializeCtrlBranch();//��ʼ��֧·��������ź�ֵ
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[3];//����ڵ��ţ�a,b,c
	int outNode[2];//����ڵ��ţ� d,q
	double inNodeValue[3];//�����ź�ֵ
	double outNodeValue[2];//����ź�ֵ
};

#endif