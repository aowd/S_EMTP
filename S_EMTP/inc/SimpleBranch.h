#ifndef SIMPLEBRANCH_H
#define SIMPLEBRANCH_H

#include "Component.h"

class SimpleBranch : public Component
{
public:
	SimpleBranch();
	~SimpleBranch(){};

	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//从节点电压数组中读支路两节点电压
	virtual void calculateBranchVoltage();//计算支路电压
	virtual void calculateBranchCurrent();//计算支路电流

	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//形成节点诺顿等效电流向量
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//形成节点导纳阵
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//保存支路电流
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//保存支路电流

	virtual void interpolate(double ratio);//给定比值，对支路的电压电流进行插值
	virtual void updateResult(int updateTypeNum);//更新开关处理过程中用于存储结果的变量
	virtual void setControlledVariableForSwitch(TVectorD ctrlStateNodeValue, int* stateNode,int nStateNode,TVectorD ctrlInputNodeValue, int* inputNode,int nInputNode){};

	// 配合平均模型的迭代校正算法使用
	// 目前仅PWMConverter和SimpleBranch中添加了这两个函数，xuyin,20121208
	virtual void storeInternalVariables(); // 存储内部变量，以便在迭代校正时恢复
	virtual void restoreInternalVariables(); // 迭代校正时恢复内部变量

public:
	double getNortonEquivalentResistance(){return nortonEquivalentResistance;};
	double getNortonEquivalentCurrent(){return nortonEquivalentCurrent;};
	double getNortonEquivalentCurrent_1(){return nortonEquivalentCurrent_1;};
	double getNortonEquivalentCurrent_2(){return nortonEquivalentCurrent_2;};
	bool isSeries();//判断支路是否有节点接地	

public:
	double nodeVoltage[2];
	double nortonEquivalentCurrent;
	double nortonEquivalentCurrent_1;
	double nortonEquivalentCurrent_2;
	double nortonEquivalentResistance;

	bool isSeriesOrNot;//存储简单元件的连接状态

	// 配合平均模型的迭代校正算法使用
	// 备份内部变量时使用
	double branchVoltage_bak;
	double branchCurrent_bak;
	double nortonEquivalentCurrent_bak;
};

#endif