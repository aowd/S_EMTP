#ifndef SINGLEPHASETRANSFORMER_H
#define SINGLEPHASETRANSFORMER_H

#include "Component.h"

class SinglePhaseTransformer:public Component {
public:
	SinglePhaseTransformer(int id,int in_nodeNumber[],double apparentPowerRating,double voltageRating[]);
	~SinglePhaseTransformer(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);
	virtual void calculateBranchVoltage();
	virtual void calculateBranchCurrent();
	virtual void calculateNortonEquivalentCurrent(double time);
	virtual void calculatePermeanceMatrix();//�������Ĵŵ�����
	virtual void calculateAdmittanceMatrix();//����ŵ�ٵ�Ч���ɾ���
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//�γɽڵ�ŵ�ٵ�Ч��������
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//�γɽڵ㵼����
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//����֧·����
	virtual void interpolate(double ratio);//������ֵ����֧·�ĵ�ѹ�������в�ֵ

private:
	double nodeVoltageArray[4];//�ڵ��ѹ	
	double nortonEquivalentCurrent[2];//ŵ�ٵ�Ч����
	double nortonEquivalentCurrent_1[2];
	double fluxArray[2];//��ͨ
	double fluxArray_1[2];

	double apparentPowerRating;//����ڹ���
	double voltageRating[2];//���ѹ
	double leakageReactance;//©�迹
	double magnetizingCurrent;//���ŵ���
	double coreAspectRatio[2];//���ĳߴ��
	double frequency;//Ƶ��

	//TMatrixD windingTurnsMatrix;//������Ȧ��������
	//TMatrixD admittanceMatrix;//ŵ�ٵ�Ч���ɾ���
	//TMatrixD permeanceMatrix;//���Ĵŵ�����

	//TMatrixD invert_MN;//�м����inv(Mss*Nss)
	//double windingTurns[2];//N1,N2,Nss={{N1,0},{0,N2}}
	//N1=V1,N2=V2
	double permeance[2];//M1,M2,Mss={{M1,M2},{M2,M1}}
	double conductance[3];//G1,G2,G3,Yss={{G1,G2},{G2,G3}}
	double invert_MN[4];//

	double branchCurrentArray[2];//֧·��ѹ
	double branchVoltageArray[2];//֧·����
	double branchVoltageArray_1[2];//��һʱ�̵�ѹ����ֵ��
	double branchCurrentArray_1[2];//��һʱ�̵�������ֵ��
};

#endif