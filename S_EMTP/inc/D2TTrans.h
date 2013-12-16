#ifndef D2TTRANS_H
#define D2TTRANS_H

#include "CtrlComponent.h"

class D2TTrans : public CtrlComponent
{
public:
	D2TTrans(int id,int nodeD,int nodeQ,int nodeA,int nodeB,int nodeC);
	~D2TTrans(){};

	virtual void initializeCtrlBranch();//��ʼ��֧·��������ź�ֵ
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[2];//����ڵ��ţ�d,q
	int outNode[3];//����ڵ��ţ� a,b,c
	double inNodeValue[2];//�����ź�ֵ
	double outNodeValue[3];//����ź�ֵ
};

#endif