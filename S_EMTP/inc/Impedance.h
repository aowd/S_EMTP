#ifndef IMPEDANCE_H
#define IMPEDANCE_H

#include "Component.h"

class Impedance:public Component
{
public:
	Impedance(int id,int fromNode,int toNode,double resistance,double inductance);
	~Impedance();
	void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流
	void readNodeVoltage(TVectorD& nodeVoltageArray);//从节点电压数组中读支路两节点电压
	void calculateBranchVoltage();//计算支路电压
	void calculateBranchCurrent();//计算支路电流
	void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	void calculateNortonEquivalentResistance(double time);//计算支路的诺顿等效电阻
	void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//形成节点诺顿等效电流向量
	void formConductanceMatrix(TMatrixD &conductanceMatrix);//形成节点导纳阵
	void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//保存支路电流
	void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//保存支路电流
	double getNortonEquivalentResistance();
	double getNortonEquivalentCurrent();
	double getNortonEquivalentCurrent_1();
	double getNortonEquivalentCurrent_2();

	void interpolate(double ratio){};//给定比值，对支路的电压电流进行插值
	void updateResult(int updateTypeNum){};//更新开关处理过程中用于存储结果的变量

	bool checkSwitch(double time){return 0;};//检测开关是否需要动作
	bool checkSwitch(int counter,TMatrixD &conductanceMatrix){return 0;};//for Inverter
	bool getState(){return 0;};//返回开关状态
	double getSwitchRatio(){return 0;};//返回开关动作点的插值比
	void switchIt(){};//变换开关状态
	void modifyConductanceMatrix(TMatrixD &conductanceMatrix){};//修正节点导纳阵
private:
	double resistance;
	double inductance;
	double *nodeVoltage;
	double nortonEquivalentResistance;
	double nortonEquivalentCurrent;
	double nortonEquivalentCurrent_1;
	double nortonEquivalentCurrent_2;
};
#endif