#ifndef INVERTER_H
#define INVERTER_H
#include "Component.h"


class Inverter:public Component{
public:
	Inverter(int id,int DCNode_first,int DCNode_last,int ACNode_first,int ACNode_last,int Inverter_number,double onValue,double offValue);
	~Inverter(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);
	virtual void calculateBranchVoltage();
	virtual void calculateBranchCurrent();
	virtual void calculateNortonEquivalentCurrent(double time);
	virtual void calculateNortonEquivalentResistance(double time);
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//保存支路电流
	virtual bool checkSwitch(int counter,TMatrixD &conductanceMatrix);
	bool switchController(double time,int Inverter_number);
	bool camparePWM(double time,int SwitchNnumber);

	void calculateIGBTVoltage();//计算各个IGBT两端的电压
	void calculateIGBTCurrent();//计算各个IGBT上通过的电流

	virtual void interpolate(double ratio);//给定比值，对支路的电压电流进行插值
	virtual void updateResult(int updateTypeNum);//更新开关处理过程中用于存储结果的变量

	void modifyConductanceMatrix(TMatrixD &conductanceMatrix,int k);//更新节点导纳矩阵
	void modifyConductanceMatrixforfirststep(TMatrixD &conductanceMatrix,int k);//更新节点导纳矩阵
	
public:
	double onValue;
	double offValue;
	double nortonEquivalentCurrent;
	double nortonEquivalentResistance[20];
	double nodeVoltage[12];
	double IGBTVoltage[20];
	double IGBTCurrent[20];
	double IGBTCurrent_1[20];
	double IGBTCurrent_2[20];
	double branchCurrentArray[20];
	double frequencyCW;
	double VmaxCW;
	double VoltageCW[20];
	double frequencyRef;
	double VmaxRef;
	double VoltageRef[20];
	int Inverter_number;
	int SwitchNnumber;
	int IGBTstate[20];

	//用于存储各个时步的状态
	//bool inverterState[10][1300];
	bool * inverterState[10];
};

#endif