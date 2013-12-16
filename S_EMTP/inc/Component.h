#ifndef COMPONENT_H
#define COMPONENT_H
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

class Component
{
public:
	Component(){isUserDef = 0;need_NEC=0;};
	virtual ~Component(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time){};//��ʼ��֧·��ѹ����
	virtual void setDetalT(double deltaT){this->deltaT = deltaT;};
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray){};//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
	virtual void calculateBranchVoltage(){};//����֧·��ѹ
	virtual void calculateBranchCurrent(){};//����֧·����
	virtual void calculateNortonEquivalentCurrent(double time){};//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time){};//����֧·��ŵ�ٵ�Ч����
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray){};//�γɽڵ�ŵ�ٵ�Ч��������
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix){};//�γɽڵ㵼����
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter){};//����֧·����
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter){};//����֧·����
	virtual void saveMachineWr(double** machineWrMatrix, int& ptr, int counter){};//������ת��
	virtual double getNortonEquivalentResistance(){return 0;};
	virtual double getNortonEquivalentCurrent(){return 0;};
	virtual double getNortonEquivalentCurrent_1(){return 0;};
	virtual double getNortonEquivalentCurrent_2(){return 0;};

	virtual void interpolate(double ratio){};//������ֵ����֧·�ĵ�ѹ�������в�ֵ
	virtual void updateResult(int updateTypeNum){};//���¿��ش�����������ڴ洢����ı���

	virtual bool checkSwitch(double time){return 0;};//��⿪���Ƿ���Ҫ����
	virtual bool checkSwitch(int counter,TMatrixD &conductanceMatrix){return 0;};//for Inverter
	virtual bool getState(){return 0;};//���ؿ���״̬
	virtual double getSwitchRatio(){return 0;};//���ؿ��ض�����Ĳ�ֵ��
	virtual void switchIt(){};//�任����״̬
	virtual void modifyConductanceMatrix(TMatrixD &conductanceMatrix){};//�����ڵ㵼����
	virtual int getSwitchMode(){return 0;};

	virtual void setControlledVariable(double* ctrlNodeValue){};
	virtual void setControlledVariableForSwitch(double* ctrlNodeValue){};

	// PWM������ר�ã�����˵���μ�PWMConverter.h,xuyin,20121208
	virtual void initializeSwitchingInstants(double** PWMInitialMatrix, int nStep, int& ptr){}; // ��ʼ������ʱ��
	virtual void ctrlNodeNumber_PWM(int* ctrlNodeNumber, int k){}; // ��¼PWM������Ӧ�Ŀ��ƽڵ���
	virtual void predictSwithcingInstants(){}; // Ԥ�⿪��ʱ��
	virtual int predictSwithcingInstants(double** PWMPredictMatrix, int nStep, int& ptr){return 0;}; // Ԥ�⿪��ʱ��
	virtual int correctSwithcingInstants(double** PWMCorrectMatrix, double tol, int nStep, int& ptr){return 0;}; // У������ʱ��
	
	// ���ƽ��ģ�͵ĵ���У���㷨ʹ��
	// Ŀǰ��PWMConverter��SimpleBranch�������������������xuyin,20121208
	virtual void storeInternalVariables(){}; // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(){}; // ����У��ʱ�ָ��ڲ�����

	// �����ʼ����xuyin, 20121226
	virtual void initializeGen(double** GenInitialMatrix, int& ptr){};
	virtual void getWindVelocityData(double** VwMatrix, int rows, int& ptr){};

	// xuyin, 20130224
	virtual void correctYne(double* tao_s){};

	//���ɱ����״̬�仯
	virtual int timedResistanceStatusChange(double time){return 0;};

public:
	int id;
	string name;
	int type;
	int nPort;
	double deltaT;
	double branchCurrent;
	double branchCurrent_1;
	double branchCurrent_2;
	double branchVoltage;
	double branchVoltage_1;
	double branchVoltage_2;
	int *nodeNumber;

	//�Ƿ��Զ���Ԫ��
	bool isUserDef;

	//�Ƿ���Ҫ����ŵ�ٵ�ֵ����
	bool need_NEC;	

	// PWM�������Ŀ��ض���ʱ��, xuyin, 20130224
	double tao_s[3];
};

#endif