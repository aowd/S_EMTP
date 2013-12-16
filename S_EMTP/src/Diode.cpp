#include "Diode.h"
#include <iostream>
using namespace std;

Diode::Diode(int id,int fromNode,int toNode,double onValue,double offValue,int state)
{
	type = 7;
	this->id = id;
	nodeNumber[0] = fromNode;
	nodeNumber[1] = toNode;
	this->onValue = onValue;
	this->offValue = offValue;
	this->state = state;
	if(state==0)
		nortonEquivalentResistance = offValue;
	else
		nortonEquivalentResistance = onValue;
	forwardVoltageDrop = 1e-5;
	//forwardVoltageDrop = 0;
	forwardBreakoverVoltage = 1e8;
	reverseWithstandVoltage = 1e8;

	isSeriesOrNot = isSeries();
}

void Diode::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent_1 = 0;
	nortonEquivalentCurrent = 0;
}

void Diode::calculateNortonEquivalentResistance(double time)
{
	if(state==0)
		nortonEquivalentResistance = offValue;
	else
		nortonEquivalentResistance = onValue;
}

void Diode::switchIt()
{
	if (state == 1)
	{
		state = 0;
		//cout<<"������"<<id<<"�ض�"<<endl;
	} 
	else
	{
		state = 1;
		//cout<<"������"<<id<<"��ͨ"<<endl;
	}
	calculateNortonEquivalentResistance(1);
}

void Diode::modifyConductanceMatrix(TMatrixD &conductanceMatrix)
{
	double temp;//����ŵ�ٵ�ֵ���ɵ�����
	if (state == 1)
	{
		temp = 1/onValue-1/offValue;
	}
	else
	{
		temp = 1/offValue-1/onValue;
	}

	int from = nodeNumber[0];
	int to = nodeNumber[1];
	if(isSeriesOrNot)
	{
		conductanceMatrix(to,to)+=temp;
		conductanceMatrix(from,from)+=temp;
		conductanceMatrix(from,to)-=temp;
		conductanceMatrix(to,from)-=temp;
	}
	else if(from==0)
		conductanceMatrix(to,to)+=temp;
	else
		conductanceMatrix(from,from)+=temp;
}

//�жϿ����Ƿ���Ҫ����
bool Diode::checkSwitch(double time)
{
	if ((state==0 && branchVoltage>forwardVoltageDrop)||(state==1 && branchCurrent<-1e-5))
//	if ((state==0 && branchVoltage>forwardVoltageDrop)||(state==1 && branchCurrent<0))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

//���ؿ��ض�����Ĳ�ֵ��
double Diode::getSwitchRatio()
{
	if (this->state == 1)
	{
		return branchCurrent_1/(branchCurrent_1-branchCurrent);
	} 
	else
	{
		return branchVoltage_1/(branchVoltage_1-branchVoltage);
	}
}

void Diode::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();
	branchCurrent = initialCurrentArray[ptr];
	ptr++;

	if (branchVoltage*branchCurrent<0)
	{
		cerr<<"����"<<id<<"��״̬�޷�ȷ����"<<endl;
	}

	//if (branchVoltage>0)
	if (branchVoltage/branchCurrent<2*onValue)
	{
		state = 1;
		calculateNortonEquivalentResistance(time);
	}
	else
	{
		state = 0;
		calculateNortonEquivalentResistance(time);
	}
}