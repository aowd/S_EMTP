#ifndef PWMCONVERTER_H
#define PWMCONVERTER_H

#include "Component.h"

#define GENERAL

class PWMConverter:public Component {

public:
	/************************************************************************/
	/*                ���캯������������                                    */
	/************************************************************************/
	PWMConverter(int firstACNode,int firstDCNode, double R, double L, double C, int firstCtrlNode, int subType);//�����ڵ�˳��: a,b,c,dc1,dc2�����ƽڵ�˳��:a,b,c,tri
	~PWMConverter();

	/************************************************************************/
	/*                ��EMTP�ӿڵĺ���                                      */
	/************************************************************************/
	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
	virtual void calculateBranchVoltage();//����֧·��ѹ
	virtual void calculateBranchCurrent();//����֧·����
	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//�γɽڵ�ŵ�ٵ�Ч��������
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//�γɽڵ㵼����
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//����֧·����

	// PWM������ƽ��ģ�����
	virtual void calculateYne(); // ����ڵ㵼�ɾ���
	virtual void initializeSwitchingInstants(double** PWMInitialMatrix, int nStep, int& ptr); // ��ʼ������ʱ��
	virtual void ctrlNodeNumber_PWM(int* ctrlNodeNumber, int k); // ��¼PWM������Ӧ�Ŀ��ƽڵ���
	
	/* Ԥ�⿪��ʱ�̣��������汾���ֱ�Ϊ�������͵���������ϵͳ*/
	/* Ŀǰʵ�ʳ����в�����һ���������ڵĿ���ʱ����ΪԤ��ֵ,��δ�������������� */
	virtual void predictSwithcingInstants(); // Ԥ�⿪��ʱ�̣��������
	virtual int predictSwithcingInstants(double** PWMPredictMatrix, int nStep, int& ptr); // Ԥ�⿪��ʱ�̣�����������ϵͳ
	
	virtual int correctSwithcingInstants(double** PWMCorrectMatrix, double tol, int nStep, int& ptr); // У������ʱ��
	
	virtual void storeInternalVariables(); // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(); // ����У��ʱ�ָ��ڲ�����

	// xuyin, 20130224
	virtual void correctYne(double* tao_s);

public:
	// ������:1,ֱ�����нӵص㣻2,ֱ����û�нӵص�
	int subType;

	// Ԫ������
	double R;
	double L;
	double C;

	// ��Ӧ���ƽڵ���
	int ctrlNodeNumber[4];

	// ŵ�ٵ�ֵ���
	double Y11[3][3];
	double Y12[3][2];
	double Y21[2][3];
	double Y22[2][2];
	double Yne[5][5];

	double Ine1[3];
	double Ine2[2];
	double Ine[5];

	// �˿ڵ�ѹ����
	double Ub1[3];
	double Ub2[2];
	double Ub[5];

	double Ib1[3];
	double Ib2[2];
	double Ib[5];
	
	// ƽ��ģ�����
	double tao_s1[3];
	double tao_s2[3];
	// double tao_s[3]; // �Ѿ��ƶ���component.h��

	double tao_s1_1[3]; // ������巨Ԥ��ʱʹ��
	double tao_s2_1[3];

	double D[3][2];

	// �����ڲ�����ʱʹ��
	double Ub1_bak[3];
	double Ub2_bak[2];
	double Ub_bak[5];

	double Ib1_bak[3];
	double Ib2_bak[2];
	double Ib_bak[5];

	//// debug, xuyin 20121223, case 279
	//int read_counter;
	//double real_tao_s[1500][3];
};

#endif