#ifndef THREEPHASETRANSFORMER_H
#define THREEPHASETRANSFORMER_H

#include "Component.h"

class ThreePhaseTransformer:public Component {
public:
	ThreePhaseTransformer(int id,int firstNode,int lastNode,double Sbase,double Vprimary,double Vsecondary,double Frequency,double Xleakage,double r_lyw,double r_ayw);
	~ThreePhaseTransformer(){};//析构函数

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//从节点电压数组中读支路两节点电压
	virtual void calculateBranchCurrent();//计算支路电流
	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void calculateNortonEquivalentResistance(double time);//计算诺顿等效导纳矩阵
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//形成节点诺顿等效电流向量
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//形成节点导纳阵
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//保存支路电流

	virtual void interpolate(double ratio);//给定比值，对支路的电压电流进行插值
	virtual void updateResult(int updateTypeNum);//更新开关处理过程中用于存储结果的变量

	virtual void storeInternalVariables(); // 存储内部变量，以便在迭代校正时恢复
	virtual void restoreInternalVariables(); // 迭代校正时恢复内部变量

private:
	int nodeNumberArray[7];//节点编号
	double Sbase;
	double Vpbase;
	double Vsbase;
	double Ibase;
	double Rbase;
	double Frequency;
	double Wbase;
	double Xleak;
	double Rlyw;
	double Alyw;
	double Imag;

	double nodeVoltage[7];
	double nodeVoltage_1[7];
	double nodeVoltage_2[7];
	double branchCurrent[7];
	double branchCurrent_1[7];
	double branchCurrent_2[7];
	double nortonEquivalentCurrent[7];
	double nortonEquivalentCurrent_1[7];
	double nortonEquivalentCurrent_2[7];
	TMatrixD nortonEquivalentConductance;

	double nodeVoltage_bak[7];
	double branchCurrent_bak[7];
	double nortonEquivalentCurrent_bak[7];
};

#endif