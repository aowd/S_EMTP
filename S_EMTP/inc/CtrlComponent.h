#ifndef CTRLCOMPONENT_H
#define CTRLCOMPONENT_H
#include <windows.h>
#ifdef WIN32
#include <w32pragma.h>
#endif
#include <Riostream.h>
#include <TMath.h>
#include <TMatrixDUtils.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompLU.h>

#include <string>

#define PI 3.141592653589793238462643383279

class CtrlComponent
{
public:
	CtrlComponent(){};
	virtual ~CtrlComponent(){};

	virtual void initializeCtrlBranch(){};//��ʼ������֧·��������ź�
	virtual void calculateCtrlEquivalentParameter(){};//����֧·�и�������
	virtual void saveInNodeValue(double* ctrlNodeValue){};//������������ź�
	virtual void saveOutNodeValue(double* ctrlNodeValue){};//������������ź�
	virtual void calculateOutputValue(double time){};//����ģ������ź�ֵ
	virtual void calculateInitOutputValue(double time){};//����ģ������ź�ֵ
	virtual int checkCalCondition(int* nodeCalMark){return 0;};
	virtual void markOutputNode(int* nodeCalMark){};
	virtual void setDeltaT(double deltaT){this->deltaT = deltaT;};

	// PWM������ƽ��ģ�����
	virtual void storeInternalVariables(){}; // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(){}; // ����У��ʱ�ָ��ڲ�����

	virtual void storeInternalVariables_pre(){}; // Ԥ��ǰ�洢
	virtual void restoreInternalVariables_pre(){}; // Ԥ���ָ�

	// xuyin, 20121110, PI��������ʼ��
	virtual void initializePICtrl(double* PICtrlInitialValue, int& ptr){};
	virtual void initializeDelay(double** DelayInitialValue, int& ptr, int nRows){};

public:
	int id;
	string name;
	int type;
	int nPort;//�˿ڸ���
	int nInPort;
	int nOutPort;
	double deltaT;
};

#endif