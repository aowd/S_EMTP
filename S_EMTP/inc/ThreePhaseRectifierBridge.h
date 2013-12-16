#ifndef THREEPHASERECTIFIERBRIDGE_H
#define THREEPHASERECTIFIERBRIDGE_H

#include "Component.h"

class ThreePhaseRectifierBridge : public Component
{
public:

	ThreePhaseRectifierBridge(int id,int *ACNode,int *DCNode,int *state,double onValue,double offValue);
	~ThreePhaseRectifierBridge(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);
	virtual void calculateBranchVoltage();
	virtual void calculateBranchCurrent();
	virtual void calculateNortonEquivalentCurrent(double time);
	virtual void calculateNortonEquivalentResistance(double time);
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);
	virtual bool checkSwitch(double time);
	virtual void interpolate(double ratio);

private:
	void numberdiode();//�������ܱ��
	void calculateDiodeVoltage();//����������������˵ĵ�ѹ
	void calculateDiodeCurrent();//���������������ͨ���ĵ���
	void calculatebranchCurrentArray();//�����ĸ�֧·�ϵĵ���
	void checkbridgemode();//�ж�������״̬��6+1��������
	void modifybridgemode();//����������״̬
	void modifydiodemode();//����������״̬
	//int getState();
	//double getratio();

private:	
	double forwardBreakoverVoltage;//���������ѹ
	double reverseWithstandVoltage;//���������ѹ
	double forwardVoltageDrop;//��ͨѹ��
	double nortonEquivalentCurrent;
	double nodeVoltage[5];
	double diodeVoltage[6];
	//double diodeVoltage_1[6];
	double diodeCurrent[6];
	//double diodeCurrent_1[6];
	double branchCurrentArray[6];
	
	double nortonEquivalentResistance[6];
	double onValue;
	double offValue;
	double ratio_1;
	double switchtime;

	int diodenodeNo[6][2];
	int diodeState[6];
//	int nodeNumber[5];//�����Ԫ���еģ�������Ҫ��ͷ�ļ��ж��������С
	int bridgeMode;
	int bridgeMode_tmp;
	int bridgeMode_1;
	int statechange[2];
	//int on_off;
};
#endif
