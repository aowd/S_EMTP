#ifndef DELAY_H
#define DELAY_H

#include "CtrlComponent.h"

class Delay : public CtrlComponent
{
public:
	Delay(int id,int inNode,int outNode,double tDelay,int nSamples);
	~Delay(){};

	virtual void initializeCtrlBranch();//��ʼ��֧·��������ź�ֵ
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual void calculateInitOutputValue(double time);
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);

	//��ʼ������
	virtual void initializeDelay(double** DelayInitialValue, int& ptr,int nRows);

	// PWM������ƽ��ģ�����
	virtual void storeInternalVariables(); // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(); // ����У��ʱ�ָ��ڲ�����

	virtual void storeInternalVariables_pre(); // Ԥ��ǰ�洢
	virtual void restoreInternalVariables_pre(); // Ԥ���ָ�

public:
	int inNode;//����ڵ���
	int outNode;//����ڵ���
	double inNodeValue;//�����ź�ֵ
	double outNodeValue;//����ź�ֵ
	double tDelay;//��ʱʱ��
	int nSamples;//������ֵ����
	int positionPointer;//����ָ��λ��
	int positionPointer_1;
	int positionPointer_2;
	int nCount;//ѭ������ָ��
	int nCount_1;
	int nCount_2;
	double* storeQueue;//������ֵ����
	double* storeQueue_1;
	double* storeQueue_2;
	double histValue;//������ʷ���ֵ
	double histValue_1;
	double histValue_2;
};

#endif