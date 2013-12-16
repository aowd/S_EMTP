#ifndef MEASURECOMPONENT_H
#define MEASURECOMPONENT_H
#include <windows.h>
#ifdef WIN32
#include <w32pragma.h>
#endif
#include <Riostream.h>
#include <TMath.h>
#include <TMatrixDUtils.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompLU.h>
#include <vector>
#include "Component.h"

class MeasureComponent
{
public:
	MeasureComponent(){};
	virtual ~MeasureComponent(){};

	virtual void initializeBranch(){measurands=0;};//��ʼ�������ź�
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches){return 0;};//��ò����ź�ֵ

public:
	int id;
	int type;
	double measurands;//�����ź�ֵ
	int ctrlNode;//��Ӧ����ϵͳ����ڵ���
};

#endif