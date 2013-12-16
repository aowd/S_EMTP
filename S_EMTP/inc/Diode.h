#ifndef DIODE_H
#define DIODE_H

#include "Breaker.h"


class Diode:public Breaker{
public:
	Diode(int id,int fromNode,int toNode,double onValue,double offValue,int state);
	~Diode(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual bool checkSwitch(double time);//��⿪�ض���
	virtual double getSwitchRatio();//���ؿ��ض�����Ĳ�ֵ��
	virtual void switchIt();//�任����״̬
	virtual void modifyConductanceMatrix(TMatrixD &conductanceMatrix);//�����ڵ㵼����
	virtual int getSwitchMode(){return 0;};


private:
	double forwardBreakoverVoltage;//���������ѹ
	double reverseWithstandVoltage;//���������ѹ
	double forwardVoltageDrop;//��ͨѹ��
};

#endif