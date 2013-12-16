#ifndef PRCOORDINATE_H
#define PRCOORDINATE_H

#include "CtrlComponent.h"

class PRCoordinate : public CtrlComponent
{
public:
	PRCoordinate(int id,int nodeA,int nodeB,int nodeC,int nodeD,int option);
	~PRCoordinate(){};

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
	int option;//ѡ�1Ϊֱ������->������ת����2Ϊ������->ֱ������
};

#endif