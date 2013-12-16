#ifndef IMPEDANCE_H
#define IMPEDANCE_H

#include "Component.h"

class Impedance:public Component
{
public:
	Impedance(int id,int fromNode,int toNode,double resistance,double inductance);
	~Impedance();
	void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����
	void readNodeVoltage(TVectorD& nodeVoltageArray);//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
	void calculateBranchVoltage();//����֧·��ѹ
	void calculateBranchCurrent();//����֧·����
	void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//�γɽڵ�ŵ�ٵ�Ч��������
	void formConductanceMatrix(TMatrixD &conductanceMatrix);//�γɽڵ㵼����
	void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//����֧·����
	void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//����֧·����
	double getNortonEquivalentResistance();
	double getNortonEquivalentCurrent();
	double getNortonEquivalentCurrent_1();
	double getNortonEquivalentCurrent_2();

	void interpolate(double ratio){};//������ֵ����֧·�ĵ�ѹ�������в�ֵ
	void updateResult(int updateTypeNum){};//���¿��ش�����������ڴ洢����ı���

	bool checkSwitch(double time){return 0;};//��⿪���Ƿ���Ҫ����
	bool checkSwitch(int counter,TMatrixD &conductanceMatrix){return 0;};//for Inverter
	bool getState(){return 0;};//���ؿ���״̬
	double getSwitchRatio(){return 0;};//���ؿ��ض�����Ĳ�ֵ��
	void switchIt(){};//�任����״̬
	void modifyConductanceMatrix(TMatrixD &conductanceMatrix){};//�����ڵ㵼����
private:
	double resistance;
	double inductance;
	double *nodeVoltage;
	double nortonEquivalentResistance;
	double nortonEquivalentCurrent;
	double nortonEquivalentCurrent_1;
	double nortonEquivalentCurrent_2;
};
#endif