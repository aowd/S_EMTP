#ifndef SAMPLEHOLD_H
#define SAMPLEHOLD_H

#include "CtrlComponent.h"

class SampleHold : public CtrlComponent
{
public:
	SampleHold(int id,int inNode, int selectNode, int outNode);
	~SampleHold(){};

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
	int SelectNode;//ѡ��ڵ���
	int outNode;//����ڵ���
	double inNodeValue;//�����ź�ֵ
	double selectNodeValue;//�����ź�ֵ
	double outNodeValue;//����ź�ֵ
	int compStatus;//0���������1������
	int compStatus_1;
	int compStatus_2;
	double histValue;//������ʷ���ֵ
	double histValue_1;
	double histValue_2;
};

#endif