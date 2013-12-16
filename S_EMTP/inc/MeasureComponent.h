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

	virtual void initializeBranch(){measurands=0;};//初始化测量信号
	virtual double getMeasurands(TVectorD nodeVoltageVec, vector<Component*>* branches){return 0;};//获得测量信号值

public:
	int id;
	int type;
	double measurands;//测量信号值
	int ctrlNode;//对应控制系统输入节点编号
};

#endif