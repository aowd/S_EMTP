#ifndef SIGMACTRLCOMPT_H
#define SIGMACTRLCOMPT_H

#include "CtrlComponent.h"

class SigmaCtrlComp : public CtrlComponent
{
public:
	SigmaCtrlComp(int id,int inNode1,int inNode2,int outNode,int option);
	~SigmaCtrlComp(){};

	virtual void initializeCtrlBranch();//��ʼ������֧·
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);
	virtual void calculateInitOutputValue(double time);

public:
	int inNode[2];//����ڵ��ż���
	int outNode;
	double inNodeValue[2];//�����ź�ֵ����
	double outNodeValue;//����ź�ֵ
	int option;//ѡ���1Ϊ�Ӻͣ�2Ϊ���3Ϊ�˻���4Ϊ����
};

#endif