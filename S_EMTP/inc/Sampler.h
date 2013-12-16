#ifndef SAMPLER_H
#define SAMPLER_H

#include "CtrlComponent.h"

class Sampler : public CtrlComponent
{
public:
	Sampler(int id,int inNode,int pulseNode, int outNode);
	~Sampler(){};

	virtual void initializeCtrlBranch();//��ʼ��֧·��������ź�ֵ
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual void calculateInitOutputValue(double time);
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);

	// PWM������ƽ��ģ�����
	virtual void storeInternalVariables(); // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(); // ����У��ʱ�ָ��ڲ�����

	virtual void storeInternalVariables_pre(); // Ԥ��ǰ�洢
	virtual void restoreInternalVariables_pre(); // Ԥ���ָ�

public:
	int inNode;//����ڵ���
	int pulseNode;
	int outNode;
	double inNodeValue;//�����ź�ֵ
	double pulseValue;//���������ź�ֵ
	double outNodeValue;//����ź�ֵ
	double histValue;
	double histValue_1; // ƽ��ģ�͵���У��ʹ��
	double histValue_2; // ƽ��ģ��Ԥ��ʹ��
};

#endif