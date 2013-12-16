#ifndef INDUCGENERATOR_H
#define INDUCGENERATOR_H

#include "Component.h"

class InducGenerator:public Component {
public:
	/************************************************************************/
	/*                构造函数与析构函数                                    */
	/************************************************************************/
	InducGenerator(int a,int b,int c);//构造函数,a、b、c为节点编号
	~InducGenerator();

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
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//保存支路电流
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//保存支路电流
	virtual void updateConductanceMatrix(TMatrixD &conductanceMatrix);//转子转速变化时更新导纳矩阵

	virtual void interpolate(double ratio);//给定比值，对支路的电压电流进行插值
	virtual void updateResult(int updateTypeNum);//更新开关处理过程中用于存储结果的变量

private:
	/************************************************************************/
	/*                发电机内部调用的函数                                  */
	/************************************************************************/
	//计算系数矩阵AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs6,Gr1...Gr5
	virtual void calculateCoefficientMatrix();
	//计算dq0坐标系下定子和转子的电压电流，并预报下一时刻定子电流
	//为下一时步计算诺顿等效电流做准备
	virtual void calculateDq0Results();
	//求解机械方程，计算转子转速
	//为下一时步求解系数矩阵、计算诺顿等值导纳做准备
	virtual void calculateOmega();
	//计算派克变换与反变换矩阵
	virtual void parkTransMatrix(TMatrixD &T,double angle);
	virtual void invParkTransMatrix(TMatrixD &T_,double angle);
public:
	/************************************************************************/
	/*                需要从外部输入的参数                                  */
	/************************************************************************/
	//用于初始化的参数
	double w;//电角速度，有名值
	double w_pu;
	//*****************
	//如何进行初始化
	//*****************

	//dq0坐标系下的电阻参数
	double Rs,Rr;
	//dq0坐标系下的电抗参数
	double Xm,Xls,Xlr;

	//定子侧电压、电流和阻抗的基值
	double Vbase,Ibase,Zbase;//其中Zbase可以计算得到

	//定转子匝数比
	double turnratio;

	//control=1,转速控制
	//control=0,转矩控制
	int control;

	//机械侧参数
	double H;
	double D;
	double Tm;//输入转矩

	/************************************************************************/
	/*				系数及系数矩阵,各矩阵的意义见相关文档                   */
	/************************************************************************/
	//
	double coff;
	double T;//T=1.0/deltaT/w*(1+coff)

	//派克变换矩阵及其逆矩阵
	TMatrixD Ps,Pr;
	TMatrixD invPs,invPr;

	//
	TMatrixD AA;//(4,4)
	TMatrixD BB;//(4,4)
	//
	TMatrixD Rss;//AA(1:2,1:2)
	TMatrixD Rsr;//AA(1:2,3:4)
	TMatrixD Rrs;//AA(3:4,1:2)
	TMatrixD Rrr;//AA(3:4,3:4)
	//
	TMatrixD Bss;//BB(1:2,1:2)
	TMatrixD Bsr;//BB(1:2,3:4)
	TMatrixD Brs;//BB(3:4,1:2)
	TMatrixD Brr;//BB(3:4,3:4)

	//处理系数矩阵，使得转换到abc坐标系时能解耦
	TMatrixD R_s;//
	TMatrixD R_ave;//
	TMatrixD R_res;//R_res=Rs-Rave;
	TMatrixD inv_Rs;//inv(R_ave)
	TMatrixD inv_Rr;//inv(Rrr)

	//
	TMatrixD Gs1;//
	TMatrixD Gs2;//
	TMatrixD Gs3;//
	TMatrixD Gs4;//
	
	//
	TMatrixD Gr1;//
	TMatrixD Gr2;//
	TMatrixD Gr3;//
	

	/************************************************************************/
	/*              各时步需要计算的物理量                                  */
	/************************************************************************/
	//以下各个物理量均为标幺值
	/*  电气侧物理量  */
	//abc坐标系下定子线圈的电压电流
	double Vsabc[3];
	double Vsabc_1[3];
	double Vsabc_2[3];
	double Isabc[3];
	double Isabc_1[3];
	double Isabc_2[3];
	//dq0坐标系下定子线圈的电压电流
	double Vsdq0[2];
	double Vsdq0_his[2];
	double Vsdq0_his2[2];
	double Isdq0[2];
	double Isdq0_his[2];//预报用
	double Isdq0_his2[2];//预报用
	double Isdq0_forcast[2];//预报值
	//abc坐标系下转子线圈的电压电流
	double Vrabc[3];
	double Vrabc_1[3];
	double Vrabc_2[3];
	double Irabc[3];
	double Irabc_1[3];
	double Irabc_2[3];
	//dq0坐标系下转子侧的电压电流
	double Vrdq0[2];
	double Vrdq0_his[2];
	double Vrdq0_his2[2];
	double Irdq0[2];
	double Irdq0_his[2];//预报用
	double Irdq0_his2[2];//预报用
	double Irdq0_forcast[2];//预报值
	/*  机械侧物理量  */
	double wr;//转子转速
	double wr_his;//wr的历史值
	double curAngle_s;//d轴领先定子a轴的角度
	double curAngle_r;//d轴领先转子a轴的角度

	/************************************************************************/
	/*              等效电路中的参数                                        */
	/************************************************************************/
	double nortonEquivalentConductance;
	double nortonEquivalentConductance_his;//等效导纳的历史值，用于更新节点导纳矩阵
	double nortonEquivalentCurrent[3];//abc坐标下的诺顿等效电流，有名值
	double nortonEquivalentCurrent_1[3];
	double nortonEquivalentCurrent_2[3];
	double nortonEquivalentCurrent_dq0[2];//dq0坐标下的诺顿等效电流，标幺值
};

#endif