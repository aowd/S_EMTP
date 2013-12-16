#ifndef TIMEDACVOLTAGESOURCE_H
#define TIMEDACVOLTAGESOURCE_H

#include "SimpleBranch.h"

class TimedAcVoltageSource:public SimpleBranch {
public:
	TimedAcVoltageSource(int id,int fromNode,int toNode,double mag,double phase,double freq,double innerResistance,double Tstart, double Tend, double dropRatio);
	~TimedAcVoltageSource(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����

private:
	double magnitude;//��ѹ��ֵ����λΪ��
	double initialAngle;//����ǣ���λΪ��
	double frequency;//Ƶ�ʣ���λΪ����
	double innerResistance;//��Դ���裬��λΪŷ
	double Tstart;//���Ͽ�ʼʱ��
	double Tend;//���Ͻ���ʱ��
	double dropRatio;//���Ϲ����е�ѹ��ԭ��ֵ֮��
};

#endif