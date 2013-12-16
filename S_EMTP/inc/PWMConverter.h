#ifndef PWMCONVERTER_H
#define PWMCONVERTER_H

#include "Component.h"

#define GENERAL

class PWMConverter:public Component {

public:
	/************************************************************************/
	/*                构造函数与析构函数                                    */
	/************************************************************************/
	PWMConverter(int firstACNode,int firstDCNode, double R, double L, double C, int firstCtrlNode, int subType);//电气节点顺序: a,b,c,dc1,dc2；控制节点顺序:a,b,c,tri
	~PWMConverter();

	/************************************************************************/
	/*                与EMTP接口的函数                                      */
	/************************************************************************/
	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//从节点电压数组中读支路两节点电压
	virtual void calculateBranchVoltage();//计算支路电压
	virtual void calculateBranchCurrent();//计算支路电流
	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//形成节点诺顿等效电流向量
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//形成节点导纳阵
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//保存支路电流

	// PWM变流器平均模型相关
	virtual void calculateYne(); // 计算节点导纳矩阵
	virtual void initializeSwitchingInstants(double** PWMInitialMatrix, int nStep, int& ptr); // 初始化开关时刻
	virtual void ctrlNodeNumber_PWM(int* ctrlNodeNumber, int k); // 记录PWM变流对应的控制节点编号
	
	/* 预测开关时刻，有两个版本，分别为线性外插和迭代求解控制系统*/
	/* 目前实际程序中采用上一个开关周期的开关时刻作为预测值,并未调用这两个函数 */
	virtual void predictSwithcingInstants(); // 预测开关时刻：线性外插
	virtual int predictSwithcingInstants(double** PWMPredictMatrix, int nStep, int& ptr); // 预测开关时刻：迭代求解控制系统
	
	virtual int correctSwithcingInstants(double** PWMCorrectMatrix, double tol, int nStep, int& ptr); // 校正开关时刻
	
	virtual void storeInternalVariables(); // 存储内部变量，以便在迭代校正时恢复
	virtual void restoreInternalVariables(); // 迭代校正时恢复内部变量

	// xuyin, 20130224
	virtual void correctYne(double* tao_s);

public:
	// 子类型:1,直流侧有接地点；2,直流侧没有接地点
	int subType;

	// 元件参数
	double R;
	double L;
	double C;

	// 对应控制节点编号
	int ctrlNodeNumber[4];

	// 诺顿等值相关
	double Y11[3][3];
	double Y12[3][2];
	double Y21[2][3];
	double Y22[2][2];
	double Yne[5][5];

	double Ine1[3];
	double Ine2[2];
	double Ine[5];

	// 端口电压电流
	double Ub1[3];
	double Ub2[2];
	double Ub[5];

	double Ib1[3];
	double Ib2[2];
	double Ib[5];
	
	// 平均模型相关
	double tao_s1[3];
	double tao_s2[3];
	// double tao_s[3]; // 已经移动到component.h中

	double tao_s1_1[3]; // 线性外插法预测时使用
	double tao_s2_1[3];

	double D[3][2];

	// 备份内部变量时使用
	double Ub1_bak[3];
	double Ub2_bak[2];
	double Ub_bak[5];

	double Ib1_bak[3];
	double Ib2_bak[2];
	double Ib_bak[5];

	//// debug, xuyin 20121223, case 279
	//int read_counter;
	//double real_tao_s[1500][3];
};

#endif