#ifndef PICTRLCOMP_H
#define PICTRLCOMP_H

#include "CtrlComponent.h"

class PICtrlComp : public CtrlComponent
{
public:
	PICtrlComp(int id,int inNode,int outNode,double pParam,double iParam);
	~PICtrlComp(){};

	virtual void initializeCtrlBranch();//��ʼ��֧·��������ź�ֵ
	virtual void calculateCtrlEquivalentParameter();//����Ԫ������
	virtual void saveInNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue);//������������ź�
	virtual void calculateOutputValue(double time);//����ģ������ź�ֵ
	virtual void calculateInitOutputValue(double time);
	virtual int checkCalCondition(int* nodeCalMark);
	virtual void markOutputNode(int* nodeCalMark);

	// PI��������ʼ��
	virtual void initializePICtrl(double* PICtrlInitialValue, int& ptr);

	// PWM������ƽ��ģ�����
	virtual void storeInternalVariables(); // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(); // ����У��ʱ�ָ��ڲ�����

	virtual void storeInternalVariables_pre(); // Ԥ��ǰ�洢
	virtual void restoreInternalVariables_pre(); // Ԥ���ָ�

public:
	int inNode;//����ڵ���
	int outNode;
	double inNodeValue;//�����ź�ֵ
	double outNodeValue;//����ź�ֵ
	double kp;//����ϵ��
	double ti;//����ϵ��  x=kp(1+1/(ti*s))u
	double c,d,c1,d1;//ģ�����
	double histCurrent;
	double histCurrent_1; // ƽ��ģ�͵���У��ʹ��
	double histCurrent_2; // ƽ��ģ��Ԥ��ʹ��
};

#endif