#ifndef NEWIGBT_H
#define NEWIGBT_H

#include "Breaker.h"


class newIGBT:public Breaker{
public:
	newIGBT(int id,int fromNode,int toNode,double onValue,double offValue,int CtrlSystemNode);
	~newIGBT(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual bool checkSwitch(double time);//��⿪�ض���
	//virtual int setControlSignal(double time);
	//virtual void readControlSignal();
	virtual double getSwitchRatio();//���ؿ��ض�����Ĳ�ֵ��
	virtual void switchIt();//�任����״̬
	virtual void modifyConductanceMatrix(TMatrixD &conductanceMatrix);//�����ڵ㵼����
	virtual int getSwitchMode();
	virtual void setControlledVariable(double* ctrlNodeValue);
	virtual void setControlledVariableForSwitch(double* ctrlNodeValue);//���ݿ����ź��趨��ֵ��ѹ

private:
	double forwardBreakoverVoltage;//���������ѹ
	double reverseWithstandVoltage;//���������ѹ
	double forwardVoltageDrop;//��ͨѹ��
	int controlSignal;
	int controlSignal_1;
	int controlSignal_2;
	int controlChangeOrNot;
	int switchMode;
	//TMatrixD controlSignalMatrix;
	int CtrlSystemNode;
};

#endif