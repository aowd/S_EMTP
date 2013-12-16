#ifndef THREEPHASETRANSFORMER_H
#define THREEPHASETRANSFORMER_H

#include "Component.h"

class ThreePhaseTransformer:public Component {
public:
	ThreePhaseTransformer(int id,int firstNode,int lastNode,double Sbase,double Vprimary,double Vsecondary,double Frequency,double Xleakage,double r_lyw,double r_ayw);
	~ThreePhaseTransformer(){};//��������

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
	virtual void calculateBranchCurrent();//����֧·����
	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����ŵ�ٵ�Ч���ɾ���
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//�γɽڵ�ŵ�ٵ�Ч��������
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//�γɽڵ㵼����
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//����֧·����

	virtual void interpolate(double ratio);//������ֵ����֧·�ĵ�ѹ�������в�ֵ
	virtual void updateResult(int updateTypeNum);//���¿��ش�����������ڴ洢����ı���

	virtual void storeInternalVariables(); // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(); // ����У��ʱ�ָ��ڲ�����

private:
	int nodeNumberArray[7];//�ڵ���
	double Sbase;
	double Vpbase;
	double Vsbase;
	double Ibase;
	double Rbase;
	double Frequency;
	double Wbase;
	double Xleak;
	double Rlyw;
	double Alyw;
	double Imag;

	double nodeVoltage[7];
	double nodeVoltage_1[7];
	double nodeVoltage_2[7];
	double branchCurrent[7];
	double branchCurrent_1[7];
	double branchCurrent_2[7];
	double nortonEquivalentCurrent[7];
	double nortonEquivalentCurrent_1[7];
	double nortonEquivalentCurrent_2[7];
	TMatrixD nortonEquivalentConductance;

	double nodeVoltage_bak[7];
	double branchCurrent_bak[7];
	double nortonEquivalentCurrent_bak[7];
};

#endif