#ifndef SINGLEPHASETRANSFORMER_H
#define SINGLEPHASETRANSFORMER_H

#include "Component.h"

class SinglePhaseTransformer:public Component {
public:
	SinglePhaseTransformer(int id,int in_nodeNumber[],double apparentPowerRating,double voltageRating[]);
	~SinglePhaseTransformer(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);
	virtual void calculateBranchVoltage();
	virtual void calculateBranchCurrent();
	virtual void calculateNortonEquivalentCurrent(double time);
	virtual void calculatePermeanceMatrix();//计算铁心磁导矩阵
	virtual void calculateAdmittanceMatrix();//计算诺顿等效导纳矩阵
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//形成节点诺顿等效电流向量
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//形成节点导纳阵
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//保存支路电流
	virtual void interpolate(double ratio);//给定比值，对支路的电压电流进行插值

private:
	double nodeVoltageArray[4];//节点电压	
	double nortonEquivalentCurrent[2];//诺顿等效电流
	double nortonEquivalentCurrent_1[2];
	double fluxArray[2];//磁通
	double fluxArray_1[2];

	double apparentPowerRating;//额定视在功率
	double voltageRating[2];//额定电压
	double leakageReactance;//漏阻抗
	double magnetizingCurrent;//励磁电流
	double coreAspectRatio[2];//铁心尺寸比
	double frequency;//频率

	//TMatrixD windingTurnsMatrix;//铁心线圈匝数矩阵
	//TMatrixD admittanceMatrix;//诺顿等效导纳矩阵
	//TMatrixD permeanceMatrix;//铁心磁导矩阵

	//TMatrixD invert_MN;//中间矩阵inv(Mss*Nss)
	//double windingTurns[2];//N1,N2,Nss={{N1,0},{0,N2}}
	//N1=V1,N2=V2
	double permeance[2];//M1,M2,Mss={{M1,M2},{M2,M1}}
	double conductance[3];//G1,G2,G3,Yss={{G1,G2},{G2,G3}}
	double invert_MN[4];//

	double branchCurrentArray[2];//支路电压
	double branchVoltageArray[2];//支路电流
	double branchVoltageArray_1[2];//上一时刻电压，插值用
	double branchCurrentArray_1[2];//上一时刻电流，插值用
};

#endif