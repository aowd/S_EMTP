#include "PhotoVoltaic.h"
#include <cmath>
#include <iostream>
using namespace std;

PhotoVoltaic::PhotoVoltaic(int id,int fromNode,int toNode,double innerResistance,int m, int n, int TctrlNode, int SctrlNode)
{
	type=26;
	this->id=id;
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	this->innerResistance = innerResistance;
	paraNumber = m;
	seriNumber = n;
	this->TctrlNode = TctrlNode;
	this->SctrlNode = SctrlNode;
	nortonEquivalentResistance = innerResistance;
	need_NEC=1;
	isSeriesOrNot = isSeries();
	//光伏阵列参数取值
	Tref = 25;// ℃
	Sref = 1000;// W/m2
	Uocr = 37;// V
	Iscr = 8.4;// A
	Umr = 29.5;// V
	Imr = 7.8;// A
	alpha = 0.0005;
	beta = 0.005;
	gamma = 0.00288;
}

void PhotoVoltaic::calculateNortonEquivalentCurrent(double time)
{
	double U = -branchVoltage/seriNumber;
	U = (U>0)? U:0;
	U = (U<37)? U:37;
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
	nortonEquivalentCurrent = Isc*(1-C1*(exp(U/C2/Uoc)-1))*paraNumber;
}

void PhotoVoltaic::calculateNortonEquivalentResistance(double time)
{
	nortonEquivalentResistance = innerResistance;
}

void PhotoVoltaic::setControlledVariable(double* ctrlNodeValue)
{
	S = ctrlNodeValue[SctrlNode-1];
	T = ctrlNodeValue[TctrlNode-1];
	calculateParameters();
}

void PhotoVoltaic::calculateParameters()
{
	double dT, dS;
	dT = T-Tref;
	dS = S/Sref-1;
	Isc=Iscr*S/Sref*(1+alpha*dT);
	Im=Imr*S/Sref*(1+alpha*dT);
	Uoc=Uocr*(1-gamma*dT)*log(exp(1.0)+beta*dS);
	Um=Umr*(1-gamma*dT)*log(exp(1.0)+beta*dS);
	C2=(Um/Uoc-1)/log(1-Im/Isc);
	C1=(1-Im/Isc)*exp(-Um/C2/Uoc);
}