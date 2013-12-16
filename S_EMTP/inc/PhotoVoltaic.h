#ifndef PHOTOVOLTAIC_H
#define PHOTOVOLTAIC_H

#include "SimpleBranch.h"

class PhotoVoltaic:public SimpleBranch {
public:
	PhotoVoltaic(int id,int formNode,int toNode,double innerResistance,int m,int n,int TctrlNode,int SctrlNode);
	~PhotoVoltaic(){};

	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual void setControlledVariable(double* ctrlNodeValue);//���ݿ����ź��趨�¶Ⱥ͹���ǿ��
	virtual void calculateParameters();//���㵱ǰ�¶������ǿ���µļ���ϵ��

private:
	double innerResistance;//����ֵ
	int paraNumber, seriNumber;//�������봮����
	int TctrlNode,SctrlNode;//�¶Ⱥ͹��տ��ƽڵ���
	double Tref, Sref;//�¶ȡ����ղο�ֵ
	double T, S;//��ǰ�¶Ⱥ͹���ֵ
	double Uocr,Iscr,Umr,Imr;//��������¶��µĿ�·��ѹ����·����������ѹ�����
	double Uoc, Isc, Um, Im;//��ǰ�������¶��µĵ�ѹ������ֵ
	double C1, C2;//����ϵ��
	double alpha, beta, gamma;//��������ϵ��

};

#endif