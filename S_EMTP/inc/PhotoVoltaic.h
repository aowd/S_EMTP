#ifndef PHOTOVOLTAIC_H
#define PHOTOVOLTAIC_H

#include "SimpleBranch.h"

class PhotoVoltaic:public SimpleBranch {
public:
	PhotoVoltaic(int id,int formNode,int toNode,double innerResistance,int m,int n,int TctrlNode,int SctrlNode);
	~PhotoVoltaic(){};

	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	virtual void setControlledVariable(double* ctrlNodeValue);//根据控制信号设定温度和光照强度
	virtual void calculateParameters();//计算当前温度与光照强度下的计算系数

private:
	double innerResistance;//内阻值
	int paraNumber, seriNumber;//并联数与串联数
	int TctrlNode,SctrlNode;//温度和光照控制节点编号
	double Tref, Sref;//温度、光照参考值
	double T, S;//当前温度和光照值
	double Uocr,Iscr,Umr,Imr;//额定光照与温度下的开路电压、短路电流、最大电压与电流
	double Uoc, Isc, Um, Im;//当前光照与温度下的电压、电流值
	double C1, C2;//计算系数
	double alpha, beta, gamma;//计算所需系数

};

#endif