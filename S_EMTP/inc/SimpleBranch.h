#ifndef SIMPLEBRANCH_H
#define SIMPLEBRANCH_H

#include "Component.h"

class SimpleBranch : public Component
{
public:
	SimpleBranch();
	~SimpleBranch(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
	virtual void calculateBranchVoltage();//����֧·��ѹ
	virtual void calculateBranchCurrent();//����֧·����

	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//�γɽڵ�ŵ�ٵ�Ч��������
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//�γɽڵ㵼����
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//����֧·����
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//����֧·����

	virtual void interpolate(double ratio);//������ֵ����֧·�ĵ�ѹ�������в�ֵ
	virtual void updateResult(int updateTypeNum);//���¿��ش�����������ڴ洢����ı���
	virtual void setControlledVariableForSwitch(TVectorD ctrlStateNodeValue, int* stateNode,int nStateNode,TVectorD ctrlInputNodeValue, int* inputNode,int nInputNode){};

	// ���ƽ��ģ�͵ĵ���У���㷨ʹ��
	// Ŀǰ��PWMConverter��SimpleBranch�������������������xuyin,20121208
	virtual void storeInternalVariables(); // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(); // ����У��ʱ�ָ��ڲ�����

public:
	double getNortonEquivalentResistance(){return nortonEquivalentResistance;};
	double getNortonEquivalentCurrent(){return nortonEquivalentCurrent;};
	double getNortonEquivalentCurrent_1(){return nortonEquivalentCurrent_1;};
	double getNortonEquivalentCurrent_2(){return nortonEquivalentCurrent_2;};
	bool isSeries();//�ж�֧·�Ƿ��нڵ�ӵ�	

public:
	double nodeVoltage[2];
	double nortonEquivalentCurrent;
	double nortonEquivalentCurrent_1;
	double nortonEquivalentCurrent_2;
	double nortonEquivalentResistance;

	bool isSeriesOrNot;//�洢��Ԫ��������״̬

	// ���ƽ��ģ�͵ĵ���У���㷨ʹ��
	// �����ڲ�����ʱʹ��
	double branchVoltage_bak;
	double branchCurrent_bak;
	double nortonEquivalentCurrent_bak;
};

#endif