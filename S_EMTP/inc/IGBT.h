#ifndef IGBT_H
#define IGBT_H

#include "Breaker.h"

class IGBT:public Breaker{
public:
	IGBT(int id,int fromNode,int toNode,double onValue,double offValue,int state,int actionTimeCase);
	~IGBT(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual bool checkSwitch(double time);//��⿪�ض���

public:
	bool switchController(double time,int SwitchNnumber);
	bool camparePWM(double time,int SwitchNnumber);

public:
	double frequencyCW;
	double VmaxCW;
	double VoltageCW;
	double frequencyRef;
	double VmaxRef;
	double VoltageRef;
	int  IGBT_number;

};

#endif